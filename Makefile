#
# Create symlink in ./test/Images to directory of raw images
# Create symlink in ./test/Processed to directory to store processed data
#
# Usage:
#	make tiff [N=1]
#   make fits [N=1]
# If mpiexec is available and N is not specified, default is to use 
# all processors.

SHELL := /bin/bash
# Raw image files ($@ expands to Makefile target, e.g., tiff, fits, jpg)
raw_dir=./test/Images/$@
# Processed .npy files
processed_dir=./test/Processed/$@

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

all:
	echo "N = $(N) PREFIX = $(PREFIX)"

tiff: $(raw_dir)
	python preProcessTIFF.py --raw_dir $(raw_dir) --processed_dir $(processed_dir)
	$(PREFIX) python fftAvg.py --processed_dir $(processed_dir)
	$(PREFIX) python specFit.py --processed_dir $(processed_dir)
	python specVis.py --processed_dir $(processed_dir)
#	python paramPlot.py --processed_dir $(processed_dir)

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