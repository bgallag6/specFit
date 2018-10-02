# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 09:49:17 2018

@author: Brendan
"""

"""
######################
# run with:
# $ mpiexec -n # python fftAvg.py    (# = number of processors)
######################
"""

"""
#######################################################
### FFT segment averaging + 3x3 Pixel Box Averaging ###
#######################################################
"""

import numpy as np
import scipy.signal
from scipy import signal
import scipy.misc
from timeit import default_timer as timer
from scipy import fftpack
import yaml
import time
import datetime
import sys
import os

havempi = True
try:
  from mpi4py import MPI
except:
  havempi = False

def fftAvg(subcube):
    
    #from scipy import fftpack
    
    pixmed = np.empty(subcube.shape[0])
    spectra_seg = np.zeros((subcube.shape[1],subcube.shape[2],len(freqs)))
    
    start_sub = timer()
    T1 = 0    
    
    for ii in range(spectra_seg.shape[0]):
        for jj in range(spectra_seg.shape[1]):        
            
            # extract timeseries + normalize by exposure time
            pixmed = subcube[:,ii,jj] / exposure     
            
            # interpolate pixel-intensity values onto specified time grid
            v_interp = np.interp(t_interp,timestamp,pixmed)  
            
            avg_array = np.zeros((len(freqs)))
            
            # trim timeseries to be integer multiple of n_segments
            v_interp = v_interp[0:len(v_interp)-rem]  
            split = np.split(v_interp, n_segments)
            
            # perform Fast Fourier Transform on each segment
            for i in range(n_segments):     
                
              sig = split[i]
              sig_fft = fftpack.fft(sig)
              #sig_fft = fftpack.rfft(sig)  # real-FFT                
              powers = np.abs(sig_fft)[pidxs]
              powers = ((powers/len(sig))**2)*(1./(sig.std()**2))*2  # normalize
              avg_array += powers
            
            avg_array /= n_segments  # average fourier power of the segments
                       
            spectra_seg[ii][jj] = np.transpose(avg_array) 
        
         
        # estimate time remaining and print to screen
        T = timer()
        T2 = T - T1
        if ii == 0:
            T_init = T - start_sub
            T_est = T_init*(spectra_seg.shape[0])  
            T_min, T_sec = divmod(T_est, 60)
            T_hr, T_min = divmod(T_min, 60)
            print("On row %i of %i, est. time remaining: %i:%.2i:%.2i" % 
                  (ii, spectra_seg.shape[0], T_hr, T_min, T_sec), flush=True)
        else:
            T_est2 = T2*(spectra_seg.shape[0]-ii)
            T_min2, T_sec2 = divmod(T_est2, 60)
            T_hr2, T_min2 = divmod(T_min2, 60)
            print("On row %i of %i, est. time remaining: %i:%.2i:%.2i" % 
                  (ii, spectra_seg.shape[0], T_hr2, T_min2, T_sec2), flush=True)
        T1 = T
             
    # print estimated and total program time to screen 
    print("Beginning est. time = %i:%.2i:%.2i" % (T_hr, T_min, T_sec), flush=True)
    T_act = timer() - start_sub
    T_min3, T_sec3 = divmod(T_act, 60)
    T_hr3, T_min3 = divmod(T_min3, 60)
    print("Actual total time = %i:%.2i:%.2i" % (T_hr3, T_min3, T_sec3), flush=True) 
                    
    return spectra_seg


"""
# Setup MPI and load datacube, time, and exposure arrays
"""  

if havempi:
  # Get_size() pulls from "-n N" specified on command line
  comm = MPI.COMM_WORLD  # set up comms
  rank = comm.Get_rank()  # Each processor gets its own "rank"
  size = comm.Get_size()
else:
  comm = None
  rank = 0
  size = 1
	
  
start = timer() 

if rank == 0:
    tStart0 = datetime.datetime.fromtimestamp(time.time())
    tStart = tStart0.strftime('%Y-%m-%d %H:%M:%S')

with open('specFit_config.yaml', 'r') as stream:
    cfg = yaml.load(stream)

directory = cfg['processed_dir']
date = cfg['date']
wavelength = cfg['wavelength']
mmap_datacube = cfg['mmap_datacube']
n_segments = cfg['num_segments']  # break data into # segments of equal length
tStep = cfg["time_step"]

if mmap_datacube == True:
    cube = np.load('%s/dataCube.npy' % directory, mmap_mode='r')
else: 
    cube = np.load('%s/dataCube.npy' % directory)

timestamp = np.load('%s/timestamps.npy' % directory)
exposure = np.load('%s/exposures.npy' % directory)

## trim top/bottom rows of cube so it divides cleanly by the # of processors
trim_top = int(np.floor((cube.shape[1] % size) / 2))
trim_bot = -int(np.ceil((cube.shape[1] % size) / 2))

if trim_top == trim_bot == 0:
    chunks = np.array_split(cube, size, axis=1)  # Split based on # processors
else:
    chunks = np.array_split(cube[:,trim_top:trim_bot], size, axis=1)

# specify which chunks should be handled by each processor
for i in range(size):
    if rank == i:
        subcube = chunks[i]

"""
if wavelength in [1600,1700]:
    time_step = 24
else:
    time_step = 12
"""

if type(tStep) == float:
    timeStep = tStep
elif type(tStep) == str:
    tDiff = list(np.diff(timestamp))
    if tStep == "min":
        timeStep = np.min(tDiff)
    elif tStep == "mode":
        timeStep = max(tDiff, key=tDiff.count)
        
# interpolate timestamps onto default-cadence time-grid
t_interp = np.linspace(0, timestamp[-1], (timestamp[-1]//timeStep)+1)  
 
# determine frequency values that FFT will evaluate   
n = len(t_interp)
rem = n % n_segments
freq_size = (n - rem) // n_segments

sample_freq = fftpack.fftfreq(freq_size, d=timeStep)
pidxs = np.where(sample_freq > 0)
freqs = sample_freq[pidxs]


## Each processor runs function on subcube, results are gathered when finished

# verify each processor received subcube with correct dimensions
ss = np.shape(subcube)
print("Processor", rank, "received an array with dimensions", ss, flush=True)

spectra_seg_part = fftAvg(subcube)

if havempi:
    spectra_seg = None 
    
    # allocate receive buffer  
    if rank == 0:
        spectra_seg = np.empty((cube.shape[1]-(trim_top-trim_bot), cube.shape[2], 
                                len(freqs)), dtype='float64')  
    
    # Gather all the results
    comm.Gather(sendbuf=spectra_seg_part, recvbuf=spectra_seg, root=0)
else:
    spectra_seg = spectra_seg_part
    
# Have one node do the last bit
if rank == 0:

    ## 3x3 Averaging
    temp = np.zeros((9,spectra_seg.shape[2]))  
    spectra_array = np.zeros((spectra_seg.shape[0]-2, spectra_seg.shape[1]-2, spectra_seg.shape[2]))
    spectra_StdDev = np.zeros((spectra_seg.shape[0]-2, spectra_seg.shape[1]-2, spectra_seg.shape[2]))
    
    # calculate 3x3 pixel-box average, start +1 & end -1 to deal with edges
    for l in range(1,spectra_seg.shape[0]-1):
        for m in range(1,spectra_seg.shape[1]-1):
            
            temp[0] = spectra_seg[l-1][m-1]
            temp[1] = spectra_seg[l-1][m]
            temp[2] = spectra_seg[l-1][m+1]
            temp[3] = spectra_seg[l][m-1]
            temp[4] = spectra_seg[l][m]
            temp[5] = spectra_seg[l][m+1]
            temp[6] = spectra_seg[l+1][m-1]
            temp[7] = spectra_seg[l+1][m]
            temp[8] = spectra_seg[l+1][m+1]
            
            spectra_array[l-1][m-1] = np.average(temp, axis=0)
            spectra_StdDev[l-1][m-1] = np.std(temp, axis=0)

    T_final = timer() - start
    T_min_final, T_sec_final = divmod(T_final, 60)
    T_hr_final, T_min_final = divmod(T_min_final, 60)
    print("Total program time = %i sec" % T_final, flush=True)
    
    print("Saving files...", flush=True)
    np.save('%s/specCube.npy' % directory, spectra_array)
    np.save('%s/specUnc.npy' % directory, spectra_StdDev)
    np.save('%s/frequencies.npy' % directory, freqs)
    
    tEnd0 = datetime.datetime.fromtimestamp(time.time())
    tEnd = tEnd0.strftime('%Y-%m-%d %H:%M:%S')
    scriptName = os.path.splitext(os.path.basename(sys.argv[0]))[0]
    
    #with open('%s_%i_region_details.txt' % (date, wavelength), 'w') as file:
    with open('log.txt', 'a+') as file:
        file.write("%s: FFT & 3x3 Averaging" % scriptName + "\n")
        file.write("----------------------------" + "\n")
        file.write("Number of time segments: %i" % n_segments + "\n")
        file.write("Program start time: %s" % tStart + "\n")
        file.write("Program end time: %s" % tEnd + "\n\n")