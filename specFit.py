# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 13:28:45 2018

@author: Brendan

Usage:
  python specFit.py --processed_dir DIR --config [specFit_config.yaml] [--mmap_spectra [True]|False]
or
  mpiexec -n N python specFit.py ...
where N = number of processors.
"""

from timeit import default_timer as timer
import numpy as np
from scipy.optimize import curve_fit as Fit
from scipy import fftpack
from scipy.stats.stats import pearsonr 
import yaml
from specModel import M1, M2         
import time
import datetime
import os   
import sys     

havempi = True
try:
  from mpi4py import MPI
except:
  havempi = False

import argparse
parser = argparse.ArgumentParser(description='specFit.py')
parser.add_argument('--processed_dir', type=str)
parser.add_argument('--mmap_spectra', type=str, default=True)
parser.add_argument('--config', type=str, default='specFit_config.yaml')
args = parser.parse_args()

processed_dir = args.processed_dir
mmap_spectra = args.mmap_spectra
config = args.config

with open(config, 'r') as stream:
    cfg = yaml.load(stream)

wavelength = cfg['wavelength']

def specFit( subcube, subcube_StdDev ):
        
  params = np.zeros((11, subcube.shape[0], subcube.shape[1]))
  
  start = timer()
  T1 = 0
   
  for l in range(subcube.shape[0]): 
    for m in range(subcube.shape[1]):
                                               
        f = freqs
        s = subcube[l][m]
        
        # use 3x3 pixel-box std.dev. or adhoc method as fitting uncertainties
        if spec_unc == 'stddev':
            ds = subcube_StdDev[l][m]     
        elif spec_unc == 'adhoc':
            ds = subcube_StdDev
                                     
        ### fit models to spectra using SciPy's Levenberg-Marquart method
        try:
            nlfit_l, nlpcov_l = Fit(M1, f, s, p0=M1_guess, bounds=(M1_low, M1_high), 
                                    sigma=ds, method='dogbox')
                  
        except RuntimeError: pass
        except ValueError: pass
    
        A, n, C = nlfit_l
        
        # first fit M2 model using 'dogbox' method          
        try:                                           
            nlfit_gp, nlpcov_gp = Fit(M2, f, s, p0=M2_guess, bounds=(M2_low, M2_high), 
                                      sigma=ds, method='dogbox', max_nfev=3000)
        
        except RuntimeError: pass
        except ValueError: pass
        
        A2, n2, C2, P2, fp2, fw2 = nlfit_gp
        
        
        # next fit M2 model using default 'trf' method
        try:            
            nlfit_gp2, nlpcov_gp2 = Fit(M2, f, s, p0=[A2, n2, C2, P2, fp2, fw2], 
                                        bounds=(M2_low, M2_high), sigma=ds, max_nfev=3000)

        except RuntimeError: pass
        except ValueError: pass
        
        A22, n22, C22, P22, fp22, fw22 = nlfit_gp2  # unpack model parameters     
                       
        # create model functions from fitted parameters
        m1_fit = M1(f, *nlfit_l)        
        m2_fit2 = M2(f, *nlfit_gp2)      
        
        #weights = subcube_StdDev[l][m]
        
        residsM1 = (s - m1_fit)
        chisqrM1 =  ((residsM1/ds)**2).sum()
        #chisqrM1 =  ((residsM1/weights)**2).sum()
        redchisqrM1 = chisqrM1 / float(f.size-3)  
        
        residsM22 = (s - m2_fit2)
        chisqrM22 = ((residsM22/ds)**2).sum()
        #chisqrM22 = ((residsM22/weights)**2).sum()
        redchisqrM22 = chisqrM22 / float(f.size-6)         
        
        f_test2 = ((chisqrM1-chisqrM22)/(6-3))/((chisqrM22)/(f.size-6))
        
        # extract the lorentzian-amplitude scaling factor
        amp_scale2 = P22 / M1(np.exp(fp22), A22, n22, C22)  
        
        
        if chisqrM1 > chisqrM22:
            rval = pearsonr(m2_fit2, s)[0]
            if C22 <= 0.:
                rollover = np.nan
            else:
                rollover = (1. / ((C22 / A22)**(-1. / n22))) / 60.
            
            # populate array with M2 parameters
            params[0][l][m] = A22
            params[1][l][m] = n22
            params[2][l][m] = C22
            params[3][l][m] = P22
            params[4][l][m] = fp22
            params[5][l][m] = fw22
            params[6][l][m] = f_test2
            params[7][l][m] = amp_scale2
            params[8][l][m] = rval
            params[9][l][m] = rollover
            params[10][l][m] = redchisqrM22
            
        else:
            rval = pearsonr(m1_fit, s)[0]
            if C <= 0.:
                rollover = np.nan
            else:
                rollover = (1. / ((C / A)**(-1. / n))) / 60.
            
            # populate array with M1 parameters
            params[0][l][m] = A
            params[1][l][m] = n
            params[2][l][m] = C
            params[3][l][m] = np.NaN
            params[4][l][m] = np.NaN
            params[5][l][m] = np.NaN
            params[6][l][m] = np.NaN
            params[7][l][m] = np.NaN
            params[8][l][m] = rval
            params[9][l][m] = rollover
            params[10][l][m] = redchisqrM1
        
        
    # estimate time remaining and print to screen
    T = timer()
    T2 = T - T1
    if l == 0:
        T_init = T - start
        T_est = T_init*(subcube.shape[0])  
        T_min, T_sec = divmod(T_est, 60)
        T_hr, T_min = divmod(T_min, 60)
        print("On row %i of %i, est. time remaining: %i:%.2i:%.2i" % 
              (l, subcube.shape[0], T_hr, T_min, T_sec), flush=True)
    else:
        T_est2 = T2*((subcube.shape[0])-l)
        T_min2, T_sec2 = divmod(T_est2, 60)
        T_hr2, T_min2 = divmod(T_min2, 60)
        print("On row %i of %i, est. time remaining: %i:%.2i:%.2i" % 
              (l, subcube.shape[0], T_hr2, T_min2, T_sec2), flush=True)
    T1 = T

  # print estimated and total program time to screen        
  print("Beginning est. time = %i:%.2i:%.2i" % (T_hr, T_min, T_sec), flush=True)
  T_act = timer() - start
  T_min3, T_sec3 = divmod(T_act, 60)
  T_hr3, T_min3 = divmod(T_min3, 60)
  print("Actual total time = %i:%.2i:%.2i" % (T_hr3, T_min3, T_sec3), flush=True) 
			
  return params
	

##############################################################################
##############################################################################

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

M1_low = cfg['M1_low']
M1_high = cfg['M1_high']
M2_low = cfg['M2_low']
M2_high = cfg['M2_high']
spec_unc = cfg['spec_unc']
M1_guess = cfg['M1_guess']
M2_guess = cfg['M2_guess']

if mmap_spectra == True:
    # load memory-mapped array as read-only
    cube = np.load('%s/specCube.npy' % processed_dir, mmap_mode='r')
    cube_StdDev = np.load('%s/specUnc.npy' % processed_dir, mmap_mode='r')
else:
    cube = np.load('%s/specCube.npy' % processed_dir)
    cube_StdDev = np.load('%s/specUnc.npy' % processed_dir)

# Split the data based on no. of processors
chunks = np.array_split(cube, size)
chunks_StdDev = np.array_split(cube_StdDev, size)

freqs = np.load('%s/frequencies.npy' % processed_dir)

# assign equal weights to all parts of curve & use as fitting uncertainties
df = np.log10(freqs[1:len(freqs)]) - np.log10(freqs[0:len(freqs)-1])
df2 = np.zeros_like(freqs)
df2[0:len(df)] = df
df2[len(df2)-1] = df2[len(df2)-2]
ds = df2

# specify which chunks should be handled by each processor
for i in range(size):
    if rank == i:
        subcube = chunks[i]
        if spec_unc == 'stddev':
            subcube_StdDev = chunks_StdDev[i]
        elif spec_unc == 'adhoc':
            subcube_StdDev = ds

# verify each processor received subcube with correct dimensions
ss = np.shape(subcube)
print("Processor", rank, "received an array with dimensions", ss, flush=True)

params_T = specFit( subcube, subcube_StdDev )

if havempi:
    # Gather all the results
    newData_p = comm.gather(params_T, root=0)

# Have one node stack the results
if rank == 0:
    if havempi:
        stack_p = np.hstack(newData_p)
    else:
        stack_p = params_T
    print(stack_p.shape, flush=True)  # Verify dimensions match input cube
 
    T_final = timer() - start
    T_min_final, T_sec_final = divmod(T_final, 60)
    T_hr_final, T_min_final = divmod(T_min_final, 60)
    print("Total program time = %i:%.2i:%.2i" % (T_hr_final, T_min_final, T_sec_final), flush=True)   
    #print("Just finished region: %s %iA" % (date, wavelength), flush=True)
  
    print("Saving files...", flush=True)
    np.save('%s/param.npy' % processed_dir, stack_p)
  
    tEnd0 = datetime.datetime.fromtimestamp(time.time())
    tEnd = tEnd0.strftime('%Y-%m-%d %H:%M:%S')
    scriptName = os.path.splitext(os.path.basename(sys.argv[0]))[0]
  
    #with open('%s_%i_region_details.txt' % (date, wavelength), 'w') as file:
    with open('log.txt', 'a+') as file:
        file.write("%s: Spectra Fitting" % scriptName + "\n")
        file.write("------------------------" + "\n")
        file.write("M1 parameter bounds (low): %s" % M1_low + "\n")
        file.write("M1 parameter bounds (high): %s" % M1_high + "\n")
        file.write("M1 initial parameter guess: %s" % M1_guess + "\n")
        file.write("M2 parameter bounds (low): %s" % M2_low + "\n")
        file.write("M2 parameter bounds (high): %s" % M2_high + "\n")
        file.write("M2 initial parameter guess: %s" % M2_guess + "\n")
        file.write("Spectra-fitting uncertainties: %s" % spec_unc + "\n")
        file.write("Program start time: %s" % tStart + "\n")
        file.write("Program end time: %s" % tEnd + "\n\n")
