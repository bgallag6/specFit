# -*- coding: utf-8 -*-
"""
Created on Fri Apr  6 15:54:30 2018

@author: Brendan
"""

"""
######################
# run with:
# $ mpiexec -n # python diff_derot_mpi.py    (# = number of processors)
######################
"""

import glob
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
import sunpy
from sunpy.map import Map
from sunpy.coordinates import frames
from sunpy.physics.solar_rotation import calculate_solar_rotate_shift
from sunpy.physics.differential_rotation import diffrot_map
from timeit import default_timer as timer
from mpi4py import MPI
import yaml
import time
import datetime
import sys
import os


def datacube(flist_chunk):
    # rebin region to desired fraction 
    def rebin(a, *args):
        shape = a.shape
        lenShape = len(shape)
        factor = np.asarray(shape)/np.asarray(args)
        evList = ['a.reshape('] + \
                 ['args[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
                 [')'] + ['.mean(%d)'%(i+1) for i in range(lenShape)]
        return eval(''.join(evList))
    
    nf1 = len(flist_chunk)
    
    exposure = np.empty((nf1))
    timestamp = np.empty((nf1))
    visAvg = np.empty((mapShape[0], mapShape[1]))
    
    # image data is int16
    dCube = np.empty((nf1, mapShape[0], mapShape[1]), dtype=np.int16)  
    
    start = timer()
    T1 = 0
    
    count = 0

    # loop through datacube and extract timeseries, timestamps, exposures
    for filename in flist_chunk:
        smap = Map(filename).submap(c3, c4)
        exposure[count] = (smap.exposure_time).value
        timestamp[count] = Time(smap.date).jd
        dmap = diffrot_map(smap, time=dt0)
        visAvg += (dmap.data / (smap.exposure_time).value)
        dCube[count] = dmap.data
        count += 1        
        
        # estimate time remaining and print to screen
        T = timer()
        T2 = T - T1
        if count == 0:
            T_init = T - start
            T_est = T_init*nf1  
            T_min, T_sec = divmod(T_est, 60)
            T_hr, T_min = divmod(T_min, 60)
            print("On row %i of %i, est. time remaining: %i:%.2i:%.2i" % 
                  (count, nf1, T_hr, T_min, T_sec), flush=True)
        else:
            T_est2 = T2*(nf1-count)
            T_min2, T_sec2 = divmod(T_est2, 60)
            T_hr2, T_min2 = divmod(T_min2, 60)
            print("On row %i of %i, est. time remaining: %i:%.2i:%.2i" % 
                  (count, nf1, T_hr2, T_min2, T_sec2), flush=True)
        T1 = T
    
    dCube_trim = dCube[:, diffLatPix:-diffLatPix, xminI:-xminF]
    visAvg = visAvg[diffLatPix:-diffLatPix, xminI:-xminF]

    np.save('%s/chunk_%i_of_%i' % (datDir, rank+1, size), dCube_trim)
    
    #return dCube_trim
    return exposure, timestamp, visAvg
  

##############################################################################
##############################################################################

comm = MPI.COMM_WORLD  # set up comms
rank = comm.Get_rank()  # Each processor gets its own "rank"

# How many processors ? (pulls from "-n 4" specified in terminal)
size = MPI.COMM_WORLD.Get_size()  

if rank == 0:
    tStart0 = datetime.datetime.fromtimestamp(time.time())
    tStart = tStart0.strftime('%Y-%m-%d %H:%M:%S')

## import program configurations from file
with open('specFit_config.yaml', 'r') as stream:
    cfg = yaml.load(stream)

imDir = cfg['fits_dir']
datDir = cfg['processed_dir']
date = cfg['date']
wavelength = cfg['wavelength']
sub_reg_coords = cfg['sub_reg_coords']

# set variables from command line
x1,x2,y1,y2 = sub_reg_coords

# create a list of all the fits files
flist = sorted(glob.glob('%s/aia*.fits' % imDir))
nf = len(flist)

# Select the middle image, to derotate around
mid_file = np.int(np.floor(nf / 2))

# make mapcube containing first and last maps & calculate derotation shifts
# if not centered at 0deg long, shift amount wont be enough -- 
# maybe only calculate latitude, use other method for longitude trim
mc_shifts = []
mapI = Map(flist[0])
mapF = Map(flist[-1])

mc_shifts.append(mapI)
mc_shifts.append(mapF)
new_mapcube1 = Map(mc_shifts, cube=True)
shifts = calculate_solar_rotate_shift(new_mapcube1, layer_index=0)

# compute longitude / latitude shift over timespan
diffLon = np.abs(np.floor((np.floor(shifts['x'][1].value)/2.)))

# calculate rotation amount in pixels
diffLonPix = diffLon / (mapI.scale)[0].value

# calculate total latitude shift in pixels, x2 since underestimated?
diffLatPix = int((np.abs((shifts['y'][1].value)) / (mapI.scale)[1].value) * 2.)
#print(diffLon, diffLonPix)
#print(diffLatPix*mapI.scale[1].value, diffLatPix)

if diffLatPix == 0:
    diffLatPix = 5  # value of zero messes with next steps
    

c3 = SkyCoord((x1-diffLon)*u.arcsec, y1*u.arcsec, frame=frames.Helioprojective)  
c4 = SkyCoord((x2+diffLon)*u.arcsec, y2*u.arcsec, frame=frames.Helioprojective) 

# get middle frame subregion & time to anchor derotation
midmap = Map(flist[mid_file]).submap(c3,c4)
dt0 = midmap.date
mapShape = midmap.data.shape

# calculate pixels to trim based off of what actual derotation trims 
# *for some reason this result is different than method above
diffMapI = diffrot_map(Map(flist[0]).submap(c3, c4), time=dt0)
diffMapF = diffrot_map(Map(flist[-1]).submap(c3, c4), time=dt0)

xminindI = np.argmin(np.fliplr(diffMapI.data), axis=1)[diffLatPix:-diffLatPix]
xminindF = np.argmin(diffMapF.data, axis=1)[diffLatPix:-diffLatPix]

xminI = mapShape[1] - np.min(xminindI)  
xminF = mapShape[1] - np.min(xminindF)


## split data and send to processors 
chunks = np.array_split(flist, size)

# specify which chunks should be handled by each processor
for i in range(size):
    if rank == i:
        subcube = chunks[i]

start = timer()

ex, t, v_avg = datacube( subcube )

# Gather all results
all_ex = comm.gather(ex, root=0)
all_t = comm.gather(t, root=0)
all_v_avg = comm.gather(v_avg, root=0)

# Have one node stack the results
if rank == 0:
    ex_arr = np.hstack(all_ex)
    tArr = np.hstack(all_t)
  
    tArr -= tArr[0]  # calculate time since first image
    tArr = np.around(tArr*86400)  # get timestamps in seconds
  
    # Calculate averaged visual image
    for j in range(len(all_v_avg)):
        if j == 0:
            vAvgArr = all_v_avg[j]
        else:
            vAvgArr += all_v_avg[j]
    vAvgArr /= nf
  
    print("Saving files...", flush=True)
    np.save('%s/exposures.npy' % datDir, ex_arr)
    np.save('%s/timestamps.npy' % datDir, tArr)
    np.save('%s/visual.npy' % datDir, vAvgArr)
  
    T_final = timer() - start
    T_min_final, T_sec_final = divmod(T_final, 60)
    T_hr_final, T_min_final = divmod(T_min_final, 60)
    print("Total program time = %i:%.2i:%.2i" % (T_hr_final, T_min_final, T_sec_final), flush=True)   
    print("Just finished region: %s %iA" % (date, wavelength), flush=True)
  
    tEnd0 = datetime.datetime.fromtimestamp(time.time())
    tEnd = tEnd0.strftime('%Y-%m-%d %H:%M:%S')
    scriptName = os.path.splitext(os.path.basename(sys.argv[0]))[0]
  
    #with open('%s_%i_region_details.txt' % (date, wavelength), 'w') as file:
    with open('log.txt', 'a+') as file:
        file.write("Region Details" + "\n")
        file.write("==============" + "\n\n")
        file.write("Date: %s" % date + "\n")
        file.write("Wavelength: %i" % wavelength + "\n\n")
        file.write("%s: Extract Data & Derotate" % scriptName + "\n")
        file.write("----------------------------------------" + "\n")
        file.write("FITS directory: %s" % imDir + "\n")
        file.write("Processed directory: %s" % datDir + "\n")
        file.write("Sub-region coordinates: (%i, %i)x, (%i, %i)y" % tuple(sub_reg_coords) + "\n")
        file.write("First image timestamp: %s" % mapI.date.strftime('%Y-%m-%d %H:%M:%S') + "\n")
        file.write("Last image timestamp: %s" % mapF.date.strftime('%Y-%m-%d %H:%M:%S') + "\n")
        file.write("Number of Images: %i" % nf + "\n")
        file.write("Program start time: %s" % tStart + "\n")
        file.write("Program end time: %s" % tEnd + "\n\n")