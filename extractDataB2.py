# -*- coding: utf-8 -*-
"""
Created on Fri Jun 22 13:31:09 2018

@author: Brendan
"""

import numpy as np
import sys
import yaml
import os

size = int(sys.argv[1])

with open('specFit_config.yaml', 'r') as stream:
    cfg = yaml.load(stream)

directory = cfg['fits_dir']
date = cfg['date']
wavelength = cfg['wavelength']
mmap_derotate = cfg['mmap_derotate']  # change to mmap_datacube
save_temp = cfg['save_temp']

cube_temp = []

# load derotated cube chunks
for i in range(size):
    temp = np.load('%s/Processed/tmp/%s/%i/chunk_%i_of_%i.npy' % (directory, date, wavelength, i+1, size))
    cube_temp.append(temp)
    
cube_final = np.vstack(cube_temp)  # stack chunks into final derotated array

del cube_temp


if mmap_derotate == "y":
    orig_shape = np.array([cube_final.shape[0], cube_final.shape[1], cube_final.shape[2]])
    
    # create memory-mapped array with similar datatype and shape to original array
    mmap_arr = np.memmap('%s/Processed/tmp/%s/%i/dataCube_mmap.npy' % (directory, date, wavelength), dtype='%s' % cube_final.dtype, mode='w+', shape=tuple(orig_shape))
    
    # write data to memory-mapped array
    mmap_arr[:] = cube_final[:]
    
    # save memory-mapped array dimensions to use when loading
    np.save('%s/Processed/tmp/%s/%i/dataCube_shape.npy' % (directory, date, wavelength), orig_shape)

    # save original array if specified
    if save_temp == "y":
        np.save('%s/Processed/tmp/%s/%i/dataCube.npy' % (directory, date, wavelength), cube_final)
        
    if save_temp == "n":
        for j in range(size):
            
            fn = '%s/Processed/tmp/%s/%i/chunk_%i_of_%i.npy' % (directory, date, wavelength, j+1, size)
            
            ## if file exists, delete it ##
            if os.path.isfile(fn):
                os.remove(fn)
            else:    ## Show an error ##
                print("Error: %s file not found" % fn)
    
    # flush memory changes to disk, then remove memory-mapped object and original array
    del mmap_arr
    del cube_final
    
else:
    np.save('%s/Processed/tmp/%s/%i/dataCube.npy' % (directory, date, wavelength), cube_final)