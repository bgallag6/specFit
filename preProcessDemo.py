# -*- coding: utf-8 -*-
"""
Usage:
  python preProcessDemo.py --processed_dir DIR --raw_dir DIR
"""
import sys
import os
from PIL import Image
import numpy as np

N = 256 # Number of images to create
prefix = 'preProcessDemo' # Image file name prefix.
Nx = 5 # Number of pixels in x 
Ny = 5 # Number of pixels in y

import argparse
parser = argparse.ArgumentParser(description='preProcessDemo.py')
parser.add_argument('--processed_dir', type=str)
parser.add_argument('--raw_dir', type=str)
args = parser.parse_args()
raw_dir = args.raw_dir
processed_dir = args.processed_dir

if not os.path.exists(raw_dir): os.makedirs(raw_dir)
if not os.path.exists(processed_dir): os.makedirs(processed_dir)

###########################################################################
# Create images
imarray = np.zeros((Nx,Ny), dtype=np.uint8)
T = float(N)
for i in range(0,N):
    arr = imarray
    t = float(i)
    tmp = np.sin(4.*2*np.pi*t/T)+np.sin(5.*2*np.pi*t/T)+np.sin(6.*2*np.pi*t/T)
    arr[:,:] = np.uint8( 255.0*(tmp/3.+1.)/2. )
    im = Image.fromarray(arr)
    im.save(os.path.join(raw_dir,prefix+'-%03d.tiff' % i))

###########################################################################
# Read images
cube = np.zeros((N,Nx,Ny), dtype=np.uint8)

timestamps = []
exposures = []
for i in range(0,N):
    # Normally we would read timestamp from Exif metadata in file
    # and convert to an integer.
    timestamps.append(i)
    
    exposures.append(1.0)

    im = Image.open(os.path.join(raw_dir,prefix+'-%03d.tiff' % i))
    cube[i,:,:] = np.asarray(im)

cube_avg = np.uint8(np.average(cube,axis=0))

###########################################################################
# Save info for specFit processing
np.save(os.path.join(processed_dir,'dataCube.npy'), cube)
np.save(os.path.join(processed_dir,'visual.npy'), cube_avg)
np.save(os.path.join(processed_dir,'timestamps.npy'), timestamps)
np.save(os.path.join(processed_dir,'exposures.npy'),exposures)

print('Wrote %d images to %s. Wrote npy files to %s.' % (N,raw_dir,processed_dir))