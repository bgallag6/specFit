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
import sunpy
from sunpy.map import Map
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames
from sunpy.physics.solar_rotation import calculate_solar_rotate_shift
from sunpy.physics.differential_rotation import diffrot_map

from timeit import default_timer as timer
from mpi4py import MPI
import yaml


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
    
    exposure = np.empty((nf1))  # exposure time
    timestamp = np.empty((nf1))  # time stamps
    vis_avg = np.zeros((midmap.data.shape[0], midmap.data.shape[1]))  # for averaged visual image
    dCube = np.empty((nf1, midmap.data.shape[0], midmap.data.shape[1]), dtype=np.int16)  # save as int16, since that is what original is
    
    start = timer()
    T1 = 0
    
    count = 0

    # loop through datacube and extract pixel data and time values
    for filename in flist_chunk:
        smap = Map(filename).submap(c3, c4)
        exposure[count] = (smap.exposure_time).value
        timestamp[count] = Time(smap.date).jd  # extract julian day time from each image
        #dmap = diffrot_map(smap, dt=-dt*u.second)
        dmap = diffrot_map(smap, time=dt0)
        vis_avg += (dmap.data / (smap.exposure_time).value)  # create normalized average visual image
        dCube[count] = dmap.data
        count += 1        
        
        # estimate time remaining and print to screen  (looks to be much better - not sure why had above?)
        T = timer()
        T2 = T - T1
        if count == 0:
            T_init = T - start
            T_est = T_init*nf1  
            T_min, T_sec = divmod(T_est, 60)
            T_hr, T_min = divmod(T_min, 60)
            print("Currently on row %i of %i, estimated time remaining: %i:%.2i:%.2i" % (count, nf1, T_hr, T_min, T_sec), flush=True)
        else:
            T_est2 = T2*(nf1-count)
            T_min2, T_sec2 = divmod(T_est2, 60)
            T_hr2, T_min2 = divmod(T_min2, 60)
            print("Currently on row %i of %i, estimated time remaining: %i:%.2i:%.2i" % (count, nf1, T_hr2, T_min2, T_sec2), flush=True)
        T1 = T
    
    dCube_trim = dCube[:,diffrot_lat_pix:-diffrot_lat_pix,xminI:-xminF]
    vis_avg = vis_avg[diffrot_lat_pix:-diffrot_lat_pix,xminI:-xminF]

    np.save('%s/Processed/tmp/%s/%i/chunk_%i_of_%i' % (directory, date, wavelength, rank+1, size), dCube_trim)
    
    #return dCube_trim
    return exposure, timestamp, vis_avg
    

comm = MPI.COMM_WORLD  # set up comms
rank = comm.Get_rank()  # Each processor gets its own "rank"
size = MPI.COMM_WORLD.Get_size()  # How many processors do we have? (pulls from "-n 4" specified in terminal execution command)


## import program configurations from file
with open('specFit_config.yaml', 'r') as stream:
    cfg = yaml.load(stream)

directory = cfg['fits_dir']
date = cfg['date']
wavelength = cfg['wavelength']
sub_reg_coords = cfg['sub_reg_coords']
coords_type = cfg['coords_type']

# set variables from command line
x1,x2,y1,y2 = sub_reg_coords

# create a list of all the fits files. This is USER-DEFINED
flist = sorted(glob.glob('%s/Images/%s/%i/aia*.fits' % (directory,date,wavelength)))
nf = len(flist)

# Select the middle image, to derotate around
mid_file = np.int(np.floor(nf / 2))

# make mapcube containing first and last maps & calculate derotation shifts
# if not centered at 0deg long, shift amount wont be enough -- maybe only calculate latitude, use other method for longitude trim
mc_shifts = []
mapI = Map(flist[0])
mapF = Map(flist[-1])

mc_shifts.append(mapI)
mc_shifts.append(mapF)
new_mapcube1 = Map(mc_shifts, cube=True)
shifts = calculate_solar_rotate_shift(new_mapcube1, layer_index=0)
diffrot_longitude = np.abs(np.floor((np.floor(shifts['x'][1].value)/2.)))  # since reflects first-last frame shift, divide by 2 for shifts around center frame
diff_long_pix = diffrot_longitude / (mapI.scale)[0].value  # calculate rotation amount in pixels
diffrot_lat_pix = int((np.abs((shifts['y'][1].value)) / (mapI.scale)[1].value) * 2.)  # calculate total latitude shift in pixels, x2 since underestimated?

if diffrot_lat_pix == 0:
    diffrot_lat_pix = 5  # value of zero messes with next steps
    

c3 = SkyCoord((x1-diffrot_longitude)*u.arcsec, y1*u.arcsec, frame=frames.Helioprojective)  
c4 = SkyCoord((x2+diffrot_longitude)*u.arcsec, y2*u.arcsec, frame=frames.Helioprojective) 

# get middle frame subregion & time to anchor derotation
midmap = Map(flist[mid_file]).submap(c3,c4)
dt0 = midmap.date

# calculate pixels to trim based off of what actual derotation trims *for some reason this result is different than method above
diff_mapI = diffrot_map(Map(flist[0]).submap(c3, c4), time=dt0)
diff_mapF = diffrot_map(Map(flist[-1]).submap(c3, c4), time=dt0)

xminindI = np.argmin(np.fliplr(diff_mapI.data),axis=1)[diffrot_lat_pix:-diffrot_lat_pix]
xminindF = np.argmin(diff_mapF.data,axis=1)[diffrot_lat_pix:-diffrot_lat_pix]

xminI = midmap.data.shape[1] - np.min(xminindI)  
xminF = midmap.data.shape[1] - np.min(xminindF)


## split data and send to processors 
chunks = np.array_split(flist, size)

# specify which chunks should be handled by each processor
for i in range(size):
    if rank == i:
        subcube = chunks[i]

start = timer()

ex, t, v_avg = datacube( subcube )  # derotate all FITS files in chunk

# Gather all results
all_ex = comm.gather(ex, root=0)
all_t = comm.gather(t, root=0)
all_v_avg = comm.gather(v_avg, root=0)

# Have one node stack the results
if rank == 0:
  ex_arr = np.hstack(all_ex)
  t_arr = np.hstack(all_t)
  
  t_arr -= t_arr[0]  # calculate time since first image
  t_arr = np.around(t_arr*86400)  # get the time value in seconds, and round to nearest whole number
  
  # Calculate averaged visual image
  for j in range(len(all_v_avg)):
      if j == 0:
          v_avg_arr = all_v_avg[j]
      else:
          v_avg_arr += all_v_avg[j]
  v_avg_arr /= nf
  
  np.save('%s/Processed/tmp/%s/%i/exposures.npy' % (directory, date, wavelength), ex_arr)
  np.save('%s/Processed/tmp/%s/%i/timestamps.npy' % (directory, date, wavelength), t_arr)
  np.save('%s/Processed/Output/%s/%i/visual.npy' % (directory, date, wavelength), v_avg_arr)
  
  T_final = timer() - start
  T_min_final, T_sec_final = divmod(T_final, 60)
  T_hr_final, T_min_final = divmod(T_min_final, 60)
  print("Total program time = %i:%.2i:%.2i" % (T_hr_final, T_min_final, T_sec_final), flush=True)   
  print("Just finished region: %s %iA" % (date, wavelength), flush=True)