# -*- coding: utf-8 -*-
"""
Usage:
  python preProcessDemo.py --processed_dir DIR --raw_dir DIR
"""
import sys
import os
from PIL import Image
import numpy as np

#N = 256 # Number of images to create
prefix = 'preProcessDemo' # Image file name prefix.
Nx = 7 # Number of pixels in x 
Ny = 7 # Number of pixels in y

import argparse
parser = argparse.ArgumentParser(description='preProcessDemo.py')
parser.add_argument('--processed_dir', type=str, default='./images/processed/demo')
parser.add_argument('--raw_dir', type=str, default='./images/raw/demo')
parser.add_argument('--Nfiles', type=int, default=256)
args = parser.parse_args()
raw_dir = args.raw_dir
processed_dir = args.processed_dir
N = args.Nfiles

if not os.path.exists(raw_dir): os.makedirs(raw_dir)
if not os.path.exists(processed_dir): os.makedirs(processed_dir)

###########################################################################
# Create images
typ = 1
imarray = np.zeros((Nx,Ny), dtype=np.uint8)
T = float(N)
for i in range(0,N):
    arr = imarray
    t = float(i)
    if typ == 1:
        # Single period
        tmp = np.sin(4*np.pi*t/T)
    elif typ == 2:
        # Power law
        f = np.arange(2,128,2)
        tmp = 0
        for j in range(0,len(f)):
            tmp = tmp + (1./float(f[j]))*np.sin(f[j]*2.*np.pi*t/T)
    else:
        # Power law + tail
        f = np.arange(2,128,2)
        tmp = 0
        for j in range(0,len(f)):
            if f[j] < 32:
                A = (1./float(f[j]))
            tmp = tmp + A*np.sin(f[j]*2.*np.pi*t/T)

    arr[:,:] = np.uint8( 255*(tmp + 1.0)/2.0 )
    im = Image.fromarray(arr)
    im.save(os.path.join(raw_dir,prefix+'-%03d.tiff' % i))

###########################################################################
# Read images
#cube = np.zeros((N,Nx,Ny), dtype=np.uint8)
cube = np.zeros((Nx,Ny,N), dtype=np.uint8)

timestamps = []
exposures = []
for i in range(0,N):
    # Normally we would read timestamp from Exif metadata in file
    # and convert to an integer.
    timestamps.append(i*2)
    
    exposures.append(1.0)

    im = Image.open(os.path.join(raw_dir,prefix+'-%03d.tiff' % i))
    #cube[i,:,:] = np.asarray(im)
    cube[:,:,i] = np.asarray(im)

#print(cube[0,0,:])

#cube_avg = np.uint8(np.average(cube,axis=0))
cube_avg = np.uint8(np.average(cube,axis=2))

###########################################################################
# Save info for specFit processing
np.save(os.path.join(processed_dir,'dataCube.npy'), cube)
np.save(os.path.join(processed_dir,'visual.npy'), cube_avg)
np.save(os.path.join(processed_dir,'timestamps.npy'), timestamps)
np.save(os.path.join(processed_dir,'exposures.npy'),exposures)

print('Wrote %d images to %s. Wrote npy files to %s.' % (N,raw_dir,processed_dir))