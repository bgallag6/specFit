# Usage:
#	make demo
# or
#   make fits
# or
#   make jpg
#
# If mpiexec is available all processors are used. To change this, specify
# the number of processors to use on the command line, e.g.
#
#   make demo N=2
#
# The directory of the raw images files can be changed using, e.g.,
#   make demo raw_dir=DIRNAME
# The directory where the processed files are placed can be changed using, e.g.,
#   make demo processed_dir=DIRNAME

SHELL := /bin/bash
# Raw image files ($@ expands to Makefile target, e.g., tiff, fits, jpg)
raw_dir=./images/raw/$@
# Processed .npy files
processed_dir=./images/processed/$@

# See if MPI is available
MPI := $(shell command -v mpiexec)
ifeq ($(MPI),)
	PREFIX=
	N=1
else
	ifeq ($(UNAME), Darwin)
		N=$(shell sysctl -n hw.physicalcpu)
	endif
	ifeq ($(UNAME), Linux)
		N=$(shell sysctl -n grep -c ^processor /proc/cpuinfo)
	endif
	PREFIX="mpiexec -n $(N)"
endif

demo:
	python preProcessDemo.py --raw_dir $(raw_dir) --processed_dir $(processed_dir)
	$(PREFIX) python fftAvg.py --processed_dir $(processed_dir)
	$(PREFIX) python specFit.py --processed_dir $(processed_dir)
	python paramPlot.py --processed_dir $(processed_dir)
	python specVis.py --processed_dir $(processed_dir)

fits:
	python preProcessFITS1.py --raw_dir $(raw_dir) --processed_dir $(processed_dir)
	python preProcessFITS2.py --processed_dir $(processed_dir)
	$(PREFIX) python fftAvg.py --processed_dir $(processed_dir)
	$(PREFIX) python specFit.py --processed_dir $(processed_dir)
	python paramPlot.py
	python specVis.py

jpg:
	python preProcessJPG1.py --raw_dir $(raw_dir) --processed_dir $(processed_dir)
	python preProcessJPG2.py --processed_dir $(processed_dir)
	$(PREFIX) python fftAvg.py --processed_dir $(processed_dir)
	$(PREFIX) python specFit.py --processed_dir $(processed_dir)
	python paramPlot.py
	python specVis.py

clean:
	rm -rf __pycache__