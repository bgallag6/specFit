---
title: 'specFit: A Python library for fitting spectral models to temporal image sequences and exploring the results'
tags:
  - Python
  - image analysis
  - power spectra
authors:
  - name: Brendan Gallagher
    orcid: 0000-0002-8353-5865
    affiliation: 1
  - name: Robert S. Weigel
    orcid: 0000-0002-9521-5228
    affiliation: 1
  - name: Karl Battams
    orcid: 0000-0002-8692-6925
    affiliation: 2
affiliations:
 - name: George Mason University
   index: 1
 - name: US Naval Research Laboratory
   index: 2
date: 10 October 2018
bibliography: paper.bib
---

# Summary

![Example figure.](https://github.com/bgallag6/specFit/blob/master/web/exampleImage1B.png)

The Atmospheric Imaging Assembly (AIA), launched in 2010 as one of three instruments on the NASA Solar Dynamics Observatory (SDO) satellite, captures one image of the solar corona every second at 4096x4096-pixel resolution [@Lemen:2012].  While this vast amount of raw image data has proven to be a valuable resource for the solar physics community, the knowledge-building possibilities of this data are currently limited by the software used to analyze it [2] (reference Karl's dissertation?). 

Animated sequences of solar images often reveal dynamic and large-scale structures that are easily visible to the eye [@Lang:2001]. Also of interest are smaller-scale periodicities in image intensity and the intensity power spectra; they can reveal information about the physical processes that are active at a given spatial location [4](get reference). Using the Fast Fourier Transform (FFT), the time series extracted from each pixel location can be transformed into the frequency domain, allowing for analysis of their spectral content.

Power spectra have been of particular use in the field of coronal seismology and in addressing the problem of coronal heating.  Work by Ireland et al looked at the average spectral properties of four regions of the corona at two different wavelengths, employing a model consisting of a power-law-plus-tail and a localized Gaussian hump [@Ireland:2015].  The authors found that each region had a characteristic power spectrum and hypothesized as to the physical drivers behind each model component.  In addition, work by Auchere et al investigated the origin of periodic signals detected in the power spectra of coronal loops, demonstrating that the signals were due to “periodic trains of pulses of random amplitudes”, rather than vibrational modes [@Auchere:2016].  With their focus on specific coronal regions or features, these studies have typically reduced the data available, either by averaging the spectra of the pixels within a region of interest, or analyzing the power content of only a targeted range of spectral frequencies.  The specFit library was developed to expand the reach of these investigations and enable routine exploration of spatially and volumetrically massive data sets.

Beyond research focused on the Sun, power spectra are investigated in extrasolar astronomical studies as well, where they provide insight into the X-ray emission variability of active galactic nuclei (AGN) [@Marshall:2015].  The break timescale in an AGN's broken-power-law power spectrum, which ``most likely corresponds with a physical timescale in [its] accretion disk'', can be used to determine whether a truncated disk is present, as well as provide support as to which physical model is the driver of the flux variability. Studies of oceanographic and atmospheric turbulence also examine power spectra to better understand rates of energy transfer [@LaCasce:2012].  Identification of a spectrum's *inertial range* can shed light on the spectral region across which energy is transferred downscale at a uniform rate, often yielding a power-law power spectrum with a $-\frac{5}{3}$ exponent, as well as the scales at which energy forcing and dissipation become important.

After a spectral model has been defined, fitting this model to an isolated spectrum is relatively straightforward.  However, when attempting to investigate any spatially large region (solar or otherwise), where each pixel represents a time series and corresponding power spectrum, the multi-step task can quickly become computationally expensive.  Through implementation of the Message Passing Interface (MPI) software (via the MPI4Py Library) in the image derotation, FFT computation, and spectra fitting steps, the ``specFit`` library allows the workload to be efficiently parallelized across an arbitrary number of computing cores.  As the image and spectra datacubes can grow to exceed the memory limitations of the MPI software, they are first memory-mapped using the Numpy package, which allows the processing of datacube sizes limited only by the system's free non-volatile memory.  

The ``specFit`` software was used in [2](Karl's dissertation) and [@Battams:2018], where it enabled the authors to investigate a 1600x1600-pixel region of the solar atmosphere over 12 hours and at 4 different wavelengths.  An initial program run-time of 100 hours per wavelength was cut to five hours by utilizing the MPI software to spread the workload across 18 computing cores.  In this study, the model parameter heatmaps were reviewed to analyze whether the various regions and structures present in the visual images had characteristic spectral properties.  Use of the ``specVis`` tool allowed the accuracy of the model fits to be checked at pixel-level resolution.

The main program is comprised of four python scripts: `preProcessFITS` `fftAvg`, `specFit`, and `paramPlot`, which are called in succession from a shell script; the number of computing cores to be used is specified using a command line switch.  

Each timeseries is linearly interpolated to satisfy the FFT requirement of uniformly-spaced data.

The ``specFit`` Python library extracts intensity time series from each pixel location in spatially co-aligned temporal sequences of 2-dimensional images, transforms each time series into power spectra via the Fast Fourier Transform, fits a parametric model to the spectra, and extracts and visualizes the best-fit parameter values as heatmaps. The included Message Passing Interface bindings allow the computationally expensive steps to be parallelized. ``specVis.py``, an included GUI-based tool, allows users to explore the spatial dependence of model parameters as heatmaps and browse pixel-level spectra and model fits.

Default parametric spectral models are included in ``specModel.py``, and instructions are provided for incorporating user-defined models. ``specFit_config.yaml`` is used to pass arguments to ``specFit.py`` such as parameter bounds.

Although this library was developed for use with extreme ultraviolet (EUV) space-based images of the solar corona [@Battams:2018], it can be used with arbitrary image sequences. Examples of both use-cases are given.
Fundamentally, all uses of this library requires a two-dimensional observation/image sequence as an input, along with a timestamp associated with each image.

Each timeseries is linearly interpolated to satisfy the FFT requirement of uniformly-spaced data.

# Acknowledgements

We acknowledge contributions from...

# References