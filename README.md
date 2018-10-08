# specFit

> A Python library for fitting spectral models to temporal image sequences and exploring the results 

![Example Image](https://github.com/bgallag6/specFit/blob/master/web/exampleImage1B.png)

The `specFit` Python library extracts intensity time series from each pixel location in spatially co-aligned temporal sequences of 2-dimensional images, transforms each time series into power spectra via the Fast Fourier Transform, fits a parametric model to the spectra, and extracts and visualizes the best-fit parameter values as heatmaps.  The included Message Passing Interface bindings allow the computationally expensive steps to be parallelized. `specVis.py`, an included GUI-based tool, allows users to explore the spatial dependence of model parameters as heatmaps and browse pixel-level spectra and model fits.

Default parametric spectral models are included in `specModel.py`, and instructions are provided for incorporating user-defined models.  The included YAML configuration file `specFit_config.yaml` is used to pass arguments to `specFit.py` such as parameter bounds.

Although this library was developed for use with extreme ultraviolet (EUV) space-based images of the solar corona [https://arxiv.org/abs/1707.02448](https://arxiv.org/abs/1707.02448), it can be used with arbitrary image sequences. Examples of both use-cases are given. 

Fundamentally, all uses of this library requires a two-dimensional observation/image sequence as an input, along with a timestamp associated with each image.

## 1. Install and Test

The dependencies of this library are the Python packages (minimum version tested): [`astropy`](https://github.com/astropy/astropy)`(2.0.3)`, [`matplotlib`](https://github.com/matplotlib/matplotlib)`(1.5.1)`, [`mpi4py`](https://github.com/mpi4py/mpi4py)`(2.0.0)`, [`numpy`](https://github.com/numpy/numpy)`(1.11.2)`, [`scipy`](https://github.com/scipy/scipy)`(0.18.1)`, [`sunpy`](https://github.com/sunpy/sunpy)`(0.8.4)`, and [`PyYAML`](https://github.com/yaml/pyyaml)`(0.1.7)`.

All of these dependencies are included in Anaconda except `sunpy`. If using Anaconda 3, the specFit library can be installed and tested using

```
git clone https://github.com/bgallag6/specFit.git
make tiff
# or
make jpg
# or
pip install sunpy && make fits
```

If `mpiexec` is installed, some of the code runs in parallel using all available processors. To force usage of only one processor, use

```
make tiff N=1
# or
make fits N=1
```

## 2. Using your own data

To use this library, a pre-processing program must create `dataCube.npy` and `timestamps.npy`. It is suggested that a user starts by modifying `preProcessDemo.py`, which generates test TIFF files and then creates the needed `.npy` files. After modifying this file, all of the processing steps can be executed using

```
make demo
```

Edit `Makefile` to see the system commands that are executed.

## 3. specVis

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
