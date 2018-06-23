# -*- coding: utf-8 -*-
"""
Created on Fri Jun 22 13:32:01 2018

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
from astropy.time import Time
import datetime
from PIL import Image
from PIL.ExifTags import TAGS

from timeit import default_timer as timer
from mpi4py import MPI
import yaml

def convert_grayscale(im_path):
  
    pic = Image.open(im_path)
    pix = np.array(pic)  # (row, column, R/G/B)

    red_image = np.array(pix[:,:,0])  # create copy of array, so wont change original
    green_image = np.array(pix[:,:,1])
    blue_image = np.array(pix[:,:,2])

    gray_image = (red_image * 0.299) + (green_image * 0.587) + (blue_image * 0.114)
    gray_image = np.around(gray_image)

    pix[:,:,0] = gray_image
    pix[:,:,1] = gray_image
    pix[:,:,2] = gray_image

#    return pix    
    return gray_image
    

def get_exif(fn):
    ret = {}
    i = Image.open(fn)
    info = i._getexif()
    for tag, value in info.items():
        decoded = TAGS.get(tag, tag)
        ret[decoded] = value
    return ret


def datacube(flist_chunk):
   
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
        dCube[count] = convert_grayscale(filename)
        exposure[count] = (smap.exposure_time).value
        vis_avg += (dmap.data / (smap.exposure_time).value)  # create normalized average visual image
        t0 = get_exif(filename)
        t1 = t0['DateTime']
        t2 = t0['DateTimeDigitized']
        timestamp[count] = datetime.datetime.strptime(t1,'%Y:%m:%d %H:%M:%S').timestamp()
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

    np.save('%s/Processed/tmp/%s/%i/chunk_%i_of_%i' % (directory, date, wavelength, rank+1, size), dCube)
    
    #return dCube_trim
    return exposure, timestamp, vis_avg
    

comm = MPI.COMM_WORLD  # set up comms
rank = comm.Get_rank()  # Each processor gets its own "rank"
size = MPI.COMM_WORLD.Get_size()  # How many processors do we have? (pulls from "-n 4" specified in terminal execution command)


## import program configurations from file
with open('specFit_config.yaml', 'r') as stream:
    cfg = yaml.load(stream)

directory = cfg['fits_dir']  # change to img_dir
date = cfg['date']
wavelength = cfg['wavelength']  # change to instr. or something
sub_reg_coords = cfg['sub_reg_coords']
coords_type = cfg['coords_type']  # get rid of

# set variables from command line
x1,x2,y1,y2 = sub_reg_coords

# create a list of all the fits files. This is USER-DEFINED
flist = sorted(glob.glob('%s/Images/%s/%i/aia*.fits' % (directory,date,wavelength)))
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
  
  np.save('%s/Processed/tmp/%s/%i/exposures.npy' % (directory, date, wavelength), ex_arr)
  np.save('%s/Processed/tmp/%s/%i/timestamps.npy' % (directory, date, wavelength), t_arr)
  np.save('%s/Processed/Output/%s/%i/visual.npy' % (directory, date, wavelength), v_avg_arr)
  
  T_final = timer() - start
  T_min_final, T_sec_final = divmod(T_final, 60)
  T_hr_final, T_min_final = divmod(T_min_final, 60)
  print("Total program time = %i:%.2i:%.2i" % (T_hr_final, T_min_final, T_sec_final), flush=True)   
  print("Just finished region: %s %iA" % (date, wavelength), flush=True)