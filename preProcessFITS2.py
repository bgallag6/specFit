# -*- coding: utf-8 -*-
"""
Created on Fri Apr  6 16:04:57 2018

@author: Brendan
"""

import numpy as np
import sys
import yaml
import os

size = int(sys.argv[1])

with open('specFit_config.yaml', 'r') as stream:
    cfg = yaml.load(stream)

directory = cfg['processed_dir']
date = cfg['date']
wavelength = cfg['wavelength']

cube_temp = []

print("Merging dataCube chunks...", flush=True)

# load derotated cube chunks
for i in range(size):
    temp = np.load('%s/chunk_%i_of_%i.npy' % (directory, i+1, size))
    cube_temp.append(temp)
    
cube_final = np.vstack(cube_temp)  # stack chunks into final derotated array

del cube_temp

print("Deleting temporary files...", flush=True)

# delete temporary chunks
for j in range(size):
    
    fn = '%s/chunk_%i_of_%i.npy' % (directory, j+1, size)
    
    ## if file exists, delete it
    if os.path.isfile(fn):
        os.remove(fn)
    else:
        print("Error: %s file not found" % fn)

print("Saving final dataCube...", flush=True)

np.save('%s/dataCube.npy' % directory, cube_final)
