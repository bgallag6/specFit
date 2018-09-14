# -*- coding: utf-8 -*-
"""
Created on Fri Aug 17 12:11:12 2018

@author: Brendan
"""





import glob
import numpy as np
from astropy.time import Time
import datetime
from PIL import Image
from PIL.ExifTags import TAGS
from timeit import default_timer as timer
from mpi4py import MPI
import yaml
#import piexif  # python 2.6+ only?


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
    dCube = np.empty((nf1, imSample.shape[0], imSample.shape[1]), dtype=np.int16)
    
    start = timer()
    T1 = 0
    
    count = 0

    # loop through datacube and extract pixel data and time/exposure values
    for filename in flist_chunk:
        dCube[count] = convertGrayscale(filename)
        timestamp[count], exposure[count] = readExif(filename)
        vis_avg += dCube[count] / exposure[count]  # create normalized average visual image
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

    np.save('%s/chunk_%i_of_%i' % (datDir, rank+1, size), dCube)
    
    #return dCube_trim
    return exposure, timestamp, vis_avg
    

comm = MPI.COMM_WORLD  # set up comms
rank = comm.Get_rank()  # Each processor gets its own "rank"
size = MPI.COMM_WORLD.Get_size()  # How many processors do we have? (pulls from "-n 4" specified in terminal execution command)

fmt = "%Y:%m:%d %H:%M:%S"

## import program configurations from file
with open('specFit_config.yaml', 'r') as stream:
    cfg = yaml.load(stream)

imDir = cfg['jpg_dir']
datDir = cfg['processed_dir']

# create a list of all the fits files. This is USER-DEFINED
flist = sorted(glob.glob('%s/*.jpg' % imDir))
nf = len(flist)

# Select the middle image, to derotate around
mid_file = np.int(np.floor(nf / 2))

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
  
  np.save('%s/exposures.npy' % datDir, ex_arr)
  np.save('%s/timestamps.npy' % datDir, t_arr)
  np.save('%s/visual.npy' % datDir, v_avg_arr)
  
  T_final = timer() - start
  T_min_final, T_sec_final = divmod(T_final, 60)
  T_hr_final, T_min_final = divmod(T_min_final, 60)
  print("Total program time = %i:%.2i:%.2i" % (T_hr_final, T_min_final, T_sec_final), flush=True)   
  #print("Just finished region: %s %iA" % (date, wavelength), flush=True)