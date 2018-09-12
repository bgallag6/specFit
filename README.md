# specFit

> A Python library for fitting and exploring spectral models to temporal image sequences



## Usage

Once inside the directory, to execute the main program, run the command:

    $ bash mainFITS.sh

The shell script will first prompt the user to specify the number of processors across which the MPI-enabled tasks will be distributed.  The individual python scipts will then be run in the following order:

1. `preProcessFITS1.py`
2. `preProcessFITS2.py`
3. `fftAvg.py`
4. `specFit.py`
5. `paramPlot.py`
6. `specVis.py`

## Dependencies

This dependencies of this library are the Python packages (minimum version tested): [`astropy`](https://github.com/astropy/astropy)`(2.0.3)`, [`matplotlib`](https://github.com/matplotlib/matplotlib)`(1.5.1)`, [`mpi4py`](https://github.com/mpi4py/mpi4py)`(2.0.0)`, [`numpy`](https://github.com/numpy/numpy)`(1.11.2)`, [`scipy`](https://github.com/scipy/scipy)`(0.18.1)`, [`sunpy`](https://github.com/sunpy/sunpy)`(0.8.4)`, [`PyYAML`](https://github.com/yaml/pyyaml)`(0.1.7)`.

## Install

To install the specFit library, run the following command:

    $ git clone https://github.com/bgallag6/specFit.git
    
