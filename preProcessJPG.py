# -*- coding: utf-8 -*-
"""
Created on Fri Aug 17 12:11:12 2018

@author: Brendan
"""

"""
######################
# run with:
# $ mpiexec -n N python preProcessJPG.py --processed_dir DIR --raw_dir DIR
# N = = number of processors
######################
"""

import glob
import numpy as np
from astropy.time import Time
import datetime
from PIL import Image
from PIL.ExifTags import TAGS
from timeit import default_timer as timer
#import piexif  # python 2.6+ only?
import os

havempi = True
try:
  from mpi4py import MPI
except:
  havempi = False

import argparse
parser = argparse.ArgumentParser(description='preProcessJPG.py')
parser.add_argument('--processed_dir', type=str)
parser.add_argument('--raw_dir', type=str)
parser.add_argument('--Nfiles', type=str, default="all")

args = parser.parse_args()

raw_dir = args.raw_dir
processed_dir = args.processed_dir
Nfiles = args.Nfiles
if not os.path.exists(raw_dir): os.makedirs(raw_dir)
if not os.path.exists(processed_dir): os.makedirs(processed_dir)


## load jpg images, convert to grayscale "matching" original FITS image   
def convertGrayscale(filename):
    im = np.array(Image.open(filename))
    redIm = np.array(im[:,:,0].astype('int16'))
    blueIm = np.array(im[:,:,1].astype('int16'))
    greenIm = np.array(im[:,:,2].astype('int16'))
    grayIm = redIm + blueIm + greenIm
    return grayIm

## extract timestamp & exposure from EXIF data
def readExif(filename):
    im = Image.open(filename)
    exif = {TAGS[k]: v for k, v in im._getexif().items()
            if k in TAGS}
    ts = datetime.datetime.strptime(exif["DateTimeOriginal"], fmt)
    timeStamp = Time(ts).jd
    expTime = float(exif["ExposureTime"][0]) / float(exif["ExposureTime"][1])
    return timeStamp, expTime


def datacube(flist_chunk):
   
    nf1 = len(flist_chunk)
    
    exposure = np.empty((nf1))  # exposure time
    timestamp = np.empty((nf1))  # time stamps
    
    imSample = convertGrayscale(flist_chunk[0])
    vis_avg = np.zeros((imSample.shape[0], imSample.shape[1]))  # for averaged visual image
    #dCube = np.empty((nf1, imSample.shape[0], imSample.shape[1]), dtype=np.int16)
    dCube = np.empty((imSample.shape[0], imSample.shape[1], nf1), dtype=np.int16)
    
    start = timer()
    T1 = 0
    
    count = 0

    # loop through datacube and extract pixel data and time/exposure values
    for filename in flist_chunk:
        #dCube[count] = convertGrayscale(filename)
        dCube[:,:,count] = convertGrayscale(filename)
        timestamp[count], exposure[count] = readExif(filename)
        #vis_avg += dCube[count] / exposure[count]  # create normalized average visual image
        vis_avg += dCube[:,:,count] / exposure[count]  # create normalized average visual image
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

    np.save('%s/chunk_%i_of_%i' % (processed_dir, rank+1, size), dCube)
    
    #return dCube_trim
    return exposure, timestamp, vis_avg

##############################################################################
if havempi:
  # Get_size() pulls from "-n N" specified on command line
  comm = MPI.COMM_WORLD   # Set up comms
  rank = comm.Get_rank()  # Each processor gets its own "rank"
  size = comm.Get_size()
else:
  comm = None
  rank = 0
  size = 1    

fmt = "%Y:%m:%d %H:%M:%S"

# create a list of all the jpg files
flist = sorted(glob.glob('%s/*.jpg' % raw_dir))

if Nfiles != "all":
    flist = flist[0:int(Nfiles)]
    
nf = len(flist)

# Select the middle image, to derotate around
mid_file = np.int(np.floor(nf / 2))

## split data and send to processors 
chunks = np.array_split(flist, size)

# specify which chunks should be handled by each processor
subcube = chunks[rank]

start = timer()

ex, t, v_avg = datacube( subcube )  # derotate all FITS files in chunk

# Gather all results
all_ex = comm.gather(ex, root=0)
all_t = comm.gather(t, root=0)
all_v_avg = comm.gather(v_avg, root=0)

# Have one node stack the results
if rank == 0:
    if havempi:
        ex_arr = np.hstack(all_ex)
        tArr = np.hstack(all_t)
    else:
        ex_arr = ex
        tArr = t
        all_v_avg = [v_avg]
  
    tArr -= tArr[0]  # calculate time since first image
    tArr = np.around(tArr*86400)  # get timestamps in seconds
  
    # Calculate averaged visual image
    for j in range(len(all_v_avg)):
        if j == 0:
            v_avg_arr = all_v_avg[j]
        else:
            v_avg_arr += all_v_avg[j]
    v_avg_arr /= nf
  
    np.save('%s/exposures.npy' % processed_dir, ex_arr)
    np.save('%s/timestamps.npy' % processed_dir, tArr)
    np.save('%s/visual.npy' % processed_dir, v_avg_arr)
  
    # Load, stack, and save dataCube chunks    
    print("Merging dataCube chunks...", flush=True)
    
    cube_temp = []
    
    # load derotated cube chunks
    for i in range(size):
        temp = np.load('%s/chunk_%i_of_%i.npy' % (processed_dir, i+1, size))
        cube_temp.append(temp)
        
    #cube_final = np.vstack(cube_temp)
    cube_final = np.dstack(cube_temp)
    
    del cube_temp
    
    # delete temporary chunks
    print("Deleting temporary files...", flush=True)
    
    for j in range(size):
        
        fn = '%s/chunk_%i_of_%i.npy' % (processed_dir, j+1, size)
        
        # if file exists, delete it
        if os.path.isfile(fn): os.remove(fn)
    
    print("Saving final dataCube...", flush=True)
    
    np.save('%s/dataCube.npy' % processed_dir, cube_final)
  
    T_final = timer() - start
    T_min_final, T_sec_final = divmod(T_final, 60)
    T_hr_final, T_min_final = divmod(T_min_final, 60)
    print("Total program time = %i:%.2i:%.2i" % (T_hr_final, T_min_final, T_sec_final), flush=True)   
    #print("Just finished region: %s %iA" % (date, wavelength), flush=True)