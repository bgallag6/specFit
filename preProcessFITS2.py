# -*- coding: utf-8 -*-
"""
Created on Fri Apr  6 16:04:57 2018

@author: Brendan
"""

import numpy as np
import sys
import yaml
import os
import glob

import argparse
parser = argparse.ArgumentParser(description='preProcessFITS2.py')
parser.add_argument('--processed_dir', type=str)
args = parser.parse_args()
processed_dir = args.processed_dir

cube_temp = []

flist = sorted(glob.glob('%s/chunk*.fits' % processed_dir))
size = len(flist)

print("Merging dataCube chunks...", flush=True)

# load derotated cube chunks
for i in range(size):
    temp = np.load('%s/chunk_%i_of_%i.npy' % (processed_dir, i+1, size))
    cube_temp.append(temp)
    
cube_final = np.vstack(cube_temp)  # stack chunks into final derotated array

del cube_temp

print("Deleting temporary files...", flush=True)

# delete temporary chunks
for j in range(size):
    
    fn = '%s/chunk_%i_of_%i.npy' % (processed_dir, j+1, size)
    
    ## if file exists, delete it
    if os.path.isfile(fn):
        os.remove(fn)
    else:
        print("Error: %s file not found" % fn)

print("Saving final dataCube...", flush=True)

np.save('%s/dataCube.npy' % processed_dir, cube_final)
