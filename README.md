# specFit

> A Python library for fitting spectral models to temporal image sequences and exploring the results 

![Example Image](https://github.com/bgallag6/specFit/blob/master/images/exampleImage1B.png)

The `specFit` Python library extracts intensity time series from each pixel location in spatially co-aligned temporal sequences of 2-dimensional images, transforms each time series into power spectra via the Fast Fourier Transform, fits a parametric model to the spectra, and extracts and visualizes the best-fit parameter values as heatmaps.  The included Message Passing Interface bindings allow the computationally expensive steps to be parallelized. `specVis.py`, an included GUI-based tool, allows users to explore the spatial dependence of model parameters as heatmaps and browse pixel-level spectra and model fits.

Example parametric spectral models are included in `specModel.py`, and instructions are provided for incorporating user-defined models.  The included YAML configuration file `specFit_config.yaml` is used to pass arguments to each script, including the directory paths for input/output data, the coordinates identifying any subregion of interest (where sub-fielding of large images is required), and flags to enable/disable memory-mapping of generated datacubes.

This library was developed for use with extreme ultraviolet (EUV) space-based images of the solar corona, but use of this library for purposes other than solar research should generally require only simple modifications to the individual scripts.  In addition, this library was developed to assist in the analysis of Flexible Image Transport System (FITS) data files; different file formats will require using a variation of the original program.  As an example, for users wishing to analyze a sequence of jpeg images, a script has been included that can extract their timestamps and exposure durations from the EXIF metadata, as well as convert from the RGB color space to grayscale to accurately capture intensity fluctuations.  

Other data types should adapt well to this software, though adjustments may be necessary depending on the user's purpose. Fundamentally, all uses of this library will require some kind of two-dimensional observation/image sequence as an input, along with their temporal (date/time) information. For those using non-FITS images, appropriate pre-processing of the observations is the responsibility of the user, and the software will assume at all times that observations are (e.g.) spatially co-aligned. 

## 1. Install and Test

### 1.1 Dependencies

The dependencies of this library are the Python packages (minimum version tested): [`astropy`](https://github.com/astropy/astropy)`(2.0.3)`, [`matplotlib`](https://github.com/matplotlib/matplotlib)`(1.5.1)`, [`mpi4py`](https://github.com/mpi4py/mpi4py)`(2.0.0)`, [`numpy`](https://github.com/numpy/numpy)`(1.11.2)`, [`scipy`](https://github.com/scipy/scipy)`(0.18.1)`, [`sunpy`](https://github.com/sunpy/sunpy)`(0.8.4)`, and [`PyYAML`](https://github.com/yaml/pyyaml)`(0.1.7)`.

To install the specFit library, run the following command:

```
git clone https://github.com/bgallag6/specFit.git
```

To test, install dependencies and then execute

```
cd specFit; bash mainFITS.sh
# or
cd specFit; bash mainJPG.sh
```

### 1.2 Usage

To execute the main program, in the top-level directory run the command:

```
bash mainFITS.sh [-n N]
```

By default, the program is run without MPI. If the number of processors `N` is given using the command line switch `-n`, the MPI-enabled tasks will be distributed across `N` processors.  The individual Python scripts will then be run in the following order:

1. `preProcessFITS1.py`
2. `preProcessFITS2.py`
3. `fftAvg.py`
4. `specFit.py`
5. `paramPlot.py`
6. `specVis.py`

## 2. specVis

To explore the `specVis` utility without having to first execute the main program, the user can edit the `specFit_config.yaml` file to read:
```
specVis_dir: "./test/validation/Processed/20120606/1600"
```
`specVis` can then be run using

```
python specVis.py
```

## 3. Contact

Brendan Gallagher <bgallag6@masonlive.gmu.edu>
