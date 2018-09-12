# specFit

> A Python library for fitting and exploring spectral models to temporal image sequences

The `specFit` Python library extracts time series from each pixel location in temporal sequences of 2-dimensional images, transforms each time series into power spectra via the Fast Fourier Transform, fits a parametetric model to the spectra, and extracts and visualizes the best-fit parameter values as heatmaps.  Included Message Passing Interface bindings allow the computationally expensive steps to be parallelized. `specVis.py`, an included GUI-based tool, allows users to explore the spatial dependence of model parameters as heatmaps and browse pixel-level spectra and model fits.

The included YAML configuration file `specFit_config.yaml` is used to pass arguments to each script.  These arguments include the directory paths from which files are loaded and to which files should be saved, the coordinates identifying the subregion of interest, and whether to memory-map the generated datacubes when loading them into the next script, among others.  The spectral models are defined in the python script, `specModel.py`.  

Use of this library for purposes other than solar research may require simple modifications to the individual scripts.  In addition, this library was developed to assist in the analysis of FITS images; different file formats will require using a variation of the original program.  For users wishing to analyze a sequence of jpeg images, a script has been included that can extract their timestamps and exposure durations from the EXIF metadata, as well as convert from the RGB color space to grayscale to accurately capture brightness fluctuations.  

Other data types should adapt well to this software, though adjustments may be necessary depending on the user's purpose. Fundamentally, all uses of this library will require some kind of two-dimensional observation/image sequence as an input, along with their temporal (date/time) information. For those using non-FITS images, appropriate pre-processing of the observations is the responsibility of the user, and the software will assume at all times that observations are (e.g.) spatially co-aligned. 

## Install

### Dependencies

The dependencies of this library are the Python packages (minimum version tested): [`astropy`](https://github.com/astropy/astropy)`(2.0.3)`, [`matplotlib`](https://github.com/matplotlib/matplotlib)`(1.5.1)`, [`mpi4py`](https://github.com/mpi4py/mpi4py)`(2.0.0)`, [`numpy`](https://github.com/numpy/numpy)`(1.11.2)`, [`scipy`](https://github.com/scipy/scipy)`(0.18.1)`, [`sunpy`](https://github.com/sunpy/sunpy)`(0.8.4)`, and [`PyYAML`](https://github.com/yaml/pyyaml)`(0.1.7)`.

To install the specFit library, run the following command:

    $ git clone https://github.com/bgallag6/specFit.git


## Usage

To execute the main program, in the top-level directory run the command:

    $ bash mainFITS.sh

The shell script will first prompt the user to specify the number of processors across which the MPI-enabled tasks will be distributed.  The individual python scipts will then be run in the following order:

1. `preProcessFITS1.py`
2. `preProcessFITS2.py`
3. `fftAvg.py`
4. `specFit.py`
5. `paramPlot.py`
6. `specVis.py`

### specVis

To explore the `specVis` utility without having to first execute the main program, the user can edit the `specFit_config.yaml` file to read:
```
specVis_dir: "./test/validation/Processed/20120606/1600"
```
`specVis` can then be run using the command

    $ python specVis.py
    
## Notes

To verify the proper functioning of the library, demo FITS and jpg images which can be processed are included in the repository under the directory `test/Images/20120606/1600`.  The script used to convert the FITS images to jpg, `fits2jpg.py`, can be found in the same location.

```
├── test
    ├── Images
    |   └── 20120606
    |       └── 1600
    |           ├── fits
    |           |   └──all_fits_images.fits
    |           └── jpg
    |               └── all_jpg_images.jpg
    ├── Processed
    |   └── 20120606
    |       └── 1600
    |           └── place_holder
    └── validation
        └── Processed
            └── 20120606
                └── 1600
                    ├── timestamps.npy
                    ├── exposures.npy
                    ├── dataCube.npy
                    ├── frequencies.npy
                    ├── specCube.npy
                    ├── specUnc.npy
                    ├── visual.npy
                    └── param.npy
```

## License
UNLICENSED
