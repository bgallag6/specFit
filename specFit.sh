#!/bin/bash

echo "The process of this program:
1) Creates datacube of raw images (solar images are derotated)
2) Extracts timeseries' of pixel intensity values with corresponding timestamps and exposure durations
3) Power-spectra are computed from extracted timeseries using the Fast Fourier Transform
4) Models are fit to the spectra and their parameters extracted"

read -p "Enter the number of processors [ex. 16]: " num

mpiexec -n $num python extractDataA1.py

python extractDataA2.py $num

mpiexec -n $num python fftAvg.py

mpiexec -n $num python specFitA.py
