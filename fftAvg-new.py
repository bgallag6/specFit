# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 09:49:17 2018

@author: Brendan

Usage:
  python fftAvg.py --processed_dir DIR [--time_step 1] [--mmap_datacube True] [--num_segments 2]
To use mpi, use
  mpiexec -n N python fftAvg.py ...
where N = number of processors
"""

desc = """
PSD segment and pixel box averaging
"""

import os
import sys
import time
import datetime
import numpy as np
from scipy import fftpack
from timeit import default_timer as timer

havempi = True
try:
    from mpi4py import MPI
except:
    havempi = False


def cfg():
    import argparse
    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--processed_dir', type=str, help='Location of input .npy files.')
    parser.add_argument('--num_segments', type=int, default=2, help='???')
    parser.add_argument('--time_step', type=str, default='mode', help='???')
    parser.add_argument('--mmap_datacube', type=bool, default=True,
                        help='mmap dataCube.npy file instead of reading into memory.')
    parser.add_argument('--box_avg_size', type=int, default=3, help='Pixel box width')
    args = parser.parse_args()

    processed_dir = args.processed_dir
    num_segments = args.num_segments
    mmap_datacube = args.mmap_datacube
    box_avg_size = args.box_avg_size
    tStep = args.time_step
    if tStep != "mode" and tStep != "min":
        try:
            tStep = float(tStep)
        except:
            raise ValueError("time_step must be a float or 'min' or 'mode'")

    if rank == 0:
        print('Using:\n' + 'processed_dir = %s,\n' % processed_dir +
              'num_segments = %s, ' % num_segments +
              'pixel-box size = %s,\n' % box_avg_size +
              'time_step = %s, ' % str(tStep) +
              'mmap_datacube = %s' % str(mmap_datacube), flush=True)

    return processed_dir, num_segments, tStep, mmap_datacube, box_avg_size


def fftAvg(subcube):
    timeseries = np.empty(subcube.shape[2])

    spectra_seg = np.zeros((subcube.shape[0], subcube.shape[1], len(freqs)))

    start_sub = timer()
    T1 = 0

    for ii in range(spectra_seg.shape[0]):
        for jj in range(spectra_seg.shape[1]):

            # extract timeseries + normalize by exposure time
            timeseries = subcube[ii, jj] / exposure

            # interpolate pixel-intensity values onto specified time grid
            v_interp = np.interp(t_interp, timestamp, timeseries)

            avg_array = np.zeros((len(freqs)))

            # trim timeseries to be integer multiple of num_segments
            v_interp = v_interp[0:len(v_interp) - rem]
            split = np.split(v_interp, num_segments)

            # perform Fast Fourier Transform on each segment
            for i in range(num_segments):
                sig = split[i]
                sig_fft = fftpack.fft(sig)
                # sig_fft = fftpack.rfft(sig)  # real-FFT
                powers = np.abs(sig_fft)[pidxs]
                if sig.std() != 0:
                    powers = ((powers / len(sig)) ** 2) * (1. / (sig.std() ** 2)) * 2  # normalize
                avg_array += powers

            avg_array /= num_segments  # average fourier power of the segments

            spectra_seg[ii][jj] = np.transpose(avg_array)

            # estimate time remaining and print to screen
        T = timer()
        T2 = T - T1
        if ii == 0:
            T_init = T - start_sub
            T_est = T_init * (spectra_seg.shape[0])
        else:
            T_est = T2 * (spectra_seg.shape[0] - ii)
        T_min, T_sec = divmod(T_est, 60)
        T_hr, T_min = divmod(T_min, 60)
        if ii == 0:
            start_time = (T_hr, T_min, T_sec)

        print("Thread %i on row %i/%i, ETR: %i:%.2i:%.2i" %
              (rank, ii, spectra_seg.shape[0], T_hr, T_min, T_sec), flush=True)
        T1 = T

    # print estimated and total program time to screen
    # TODO: What is point of this next print statement?
    print("Beginning est. time = %i:%.2i:%.2i" % start_time, flush=True)
    T_act = timer() - start_sub
    T_min, T_sec = divmod(T_act, 60)
    T_hr, T_min = divmod(T_min, 60)
    print("Actual total time = %i:%.2i:%.2i" % (T_hr, T_min, T_sec), flush=True)

    return spectra_seg


"""
Set-up MPI and load datacube, time, and exposure arrays from .npy files
"""

if havempi:
    # TODO: What if mpi4py is found but user did not launch using
    #   mpiexec -n N python ...
    # ?
    # Get_size() pulls from "-n N" specified on command line
    comm = MPI.COMM_WORLD  # Set up comms
    rank = comm.Get_rank()  # Each processor gets its own "rank"
    size = comm.Get_size()
else:
    comm = None
    rank = 0
    size = 1

(processed_dir, num_segments, tStep, mmap_datacube, box_avg_size) = cfg()

start = timer()

if rank == 0:
    tStart0 = datetime.datetime.fromtimestamp(time.time())
    tStart = tStart0.strftime('%Y-%m-%d %H:%M:%S')

try:
    if mmap_datacube:
        cube = np.load('%s/dataCube.npy' % processed_dir, mmap_mode='r')
    else:
        cube = np.load('%s/dataCube.npy' % processed_dir)

    timestamp = np.load('%s/timestamps.npy' % processed_dir)
    exposure = np.load('%s/exposures.npy' % processed_dir)
    vis0 = np.load('%s/visual.npy' % processed_dir)
except FileNotFoundError:
    if rank == 0:
        sys.exit("One or more required input files not found. Exiting.")
    else:
        sys.exit(1)

# TODO: Warn if any trimming performed.
# trim top/bottom rows of cube so it divides cleanly by the # of processors
trim_top = int(np.floor((cube.shape[0] % size) / 2))
trim_bot = -int(np.ceil((cube.shape[0] % size) / 2))

# trim region border to account for pixel-box averaging
box_trim = ((box_avg_size - 1) // 2)

if rank == 0:
    print("Reading time-averaged image file %s/visual-trimmed.npy" % processed_dir)
    vis = vis0[trim_top + box_trim:vis0.shape[0] + trim_bot - box_trim, box_trim:vis0.shape[1] - box_trim]
    # TODO: This will over-write previously computed visual.npy.
    #       Give it a different name, e.g., visual-trimmed.npy?
    np.save('%s/visual.npy' % processed_dir, vis)
    print("Saving trimmed time-averaged image file %s/visual-trimmed.npy" % processed_dir)

# TODO: Remove this line or add comment for why it is commented out?
# chunks = np.array_split(cube[trim_top:cube.shape[0]+trim_bot], size, axis=0)

# TODO: Remove this line or add comment for why it is commented out?
# subcube = chunks[rank]
subcube = np.array_split(cube[trim_top:cube.shape[0] + trim_bot], size, axis=0)[rank]

if type(tStep) == float:
    timeStep = tStep
elif type(tStep) == str:
    tDiff = list(np.diff(timestamp))
    if tStep == "min":
        timeStep = np.min(tDiff)
        if rank == 0:
            print("Minimum time difference between timestamps in timestamps.npy is %d" % timeStep)
    elif tStep == "mode":
        timeStep = max(tDiff, key=tDiff.count)
        if rank == 0:
            print("Mode of time difference between timestamps in timestamps.npy is %.1f" % timeStep)

# interpolate timestamps onto default-cadence time-grid
t_interp = np.linspace(0, timestamp[-1], (timestamp[-1] // timeStep) + 1)

# determine frequency values that FFT will evaluate
n = len(t_interp)
rem = n % num_segments
freq_size = (n - rem) // num_segments

sample_freq = fftpack.fftfreq(freq_size, d=timeStep)

pidxs = np.where(sample_freq > 0)

if freq_size % 2 == 0:  # Even time series length. Keep f = -0.5 value.
    pidxs = np.append(pidxs, [pidxs[0][-1] + 1])

freqs = sample_freq[pidxs]
if freq_size % 2 == 0:
    freqs[-1] = -freqs[-1]

# Each processor computes and averages `num_segments` spectra for each pixel
# in subcube; results are combined when all processing is complete.

# Report when a processor receives a subcube
ss = np.shape(subcube)
print("Processor", rank, "received an array with dimensions", ss, flush=True)

spectra_seg_part = fftAvg(subcube)

if havempi:
    spectra_seg = None

    # allocate receive buffer
    if rank == 0:
        spectra_seg = np.empty((cube.shape[0] - (trim_top - trim_bot), cube.shape[1],
                                len(freqs)), dtype='float64')

    # Gather all the results
    comm.Gather(sendbuf=spectra_seg_part, recvbuf=spectra_seg, root=0)
else:
    spectra_seg = spectra_seg_part

# Have one node do pixel-box averaging, save output files, and write log file.
if rank == 0:

    ## Pixel-box averaging
    temp = np.zeros((box_avg_size ** 2, spectra_seg.shape[2]))
    spectra_array = np.zeros(
        (spectra_seg.shape[0] - (box_trim * 2), spectra_seg.shape[1] - (box_trim * 2), spectra_seg.shape[2]))
    if box_avg_size > 1:
        spectra_StdDev = np.zeros(
            (spectra_seg.shape[0] - (box_trim * 2), spectra_seg.shape[1] - (box_trim * 2), spectra_seg.shape[2]))

    # calculate pixel-box average, modified range to deal with edges
    for l in range(box_trim, spectra_seg.shape[0] - box_trim):
        for m in range(box_trim, spectra_seg.shape[1] - box_trim):

            count = 0

            for l1 in range(-box_trim, box_trim + 1, 1):
                for m2 in range(-box_trim, box_trim + 1, 1):
                    temp[count] = spectra_seg[l + l1][m + m2]
                    count += 1

            spectra_array[l - box_trim][m - box_trim] = np.average(temp, axis=0)
            if box_avg_size > 1:
                spectra_StdDev[l - box_trim][m - box_trim] = np.std(temp, axis=0)

    T_final = timer() - start
    T_min_final, T_sec_final = divmod(T_final, 60)
    T_hr_final, T_min_final = divmod(T_min_final, 60)
    print("Total calculation time = %i sec" % T_final, flush=True)

    print("Saving %s/specCube.npy" % processed_dir, flush=True)
    np.save('%s/specCube.npy' % processed_dir, spectra_array)

    print("Saving %s/frequencies.npy" % processed_dir, flush=True)
    np.save('%s/frequencies.npy' % processed_dir, freqs)

    if box_avg_size > 1:
        print("Saving %s/specUnc.npy" % processed_dir, flush=True)
        np.save('%s/specUnc.npy' % processed_dir, spectra_StdDev)

    tEnd0 = datetime.datetime.fromtimestamp(time.time())
    tEnd = tEnd0.strftime('%Y-%m-%d %H:%M:%S')
    scriptName = os.path.splitext(os.path.basename(sys.argv[0]))[0]

    # TODO: Save log file to same directory as output files
    # TODO: Why is this appending? Would think one run = one log file.
    with open('log.txt', 'a+') as f:
        f.write("%s: FFT & Pixel-Box Averaging" % scriptName + "\n")
        f.write("----------------------------" + "\n")
        f.write("Number of time segments: %i" % num_segments + "\n")
        f.write("Pixel-box size: %i" % box_avg_size + "\n")
        f.write("Program start time: %s" % tStart + "\n")
        f.write("Program end time: %s" % tEnd + "\n\n")