# Image Degradation
MATLAB code for reducing spatial resolution of histopathology images. Image datasets used:
* BACH (data [here](https://iciar2018-challenge.grand-challenge.org/), paper [here](https://www.sciencedirect.com/science/article/pii/S1361841518307941))
* BreaKHis (data [here](https://web.inf.ufpr.br/vri/databases/breast-cancer-histopathological-database-breakhis/), paper [here](https://ieeexplore.ieee.org/abstract/document/7312934/))
* PatchCamelyon (data [here](https://www.kaggle.com/c/histopathologic-cancer-detection), paper [here](https://link.springer.com/chapter/10.1007/978-3-030-00934-2_24)) 

## Functions
* `dataset_parameters.m`: returns struct called `params` containing previously determined parameters of imaging system used to capture publicly available datasets. Takes a single required argument:
    * `dataset` - (required) string, one of 'BreaKHis4X', 'BreaKHis10X', 'BreaKHis20X', 'BreaKHis40X', 'BACH', 'CAMELYON16-UMCU', 'CAMELYON16-RUMC', 'PCam'
    * `params` struct contains `NA_H` (original NA of system), `M_obj` (magnification of objective), `M_relay` (magnification of any relay used, 1 otherwise, `d_pixel` (physical size of camera pixel in metres), `f` (focal length of objective), `img_ext` (expected extension of image files in dataset), `scale_factor` (integer number by which to upscale images if necessary, e.g. for PCam)
* `NA_estimation.m`: code demonstrating process used to estimate NA_H for all datasets. While this calls `dataset_parameters`, it does not use the stored value of `NA_H` for that dataset for computation. It is overlaid on the resulting plot to illustrate how the value was obtained for each dataset. Takes the following arguments:
    * `dataset` - (required) string, one of 'BreaKHis4X', 'BreaKHis10X', 'BreaKHis20X', 'BreaKHis40X', 'BACH', 'CAMELYON16-UMCU', 'CAMELYON16-RUMC', 'PCam'
    * `cutoff_min` - (optional, default 10) integer, divisor of camera sampling frequency `f_s` to set minimum cutoff frequency for `NA_H` estimation
    * `show_estimation` - (optional, default `false`) boolean flag of whether to display images during NA estimation (for debugging)
* `NA_degradation.m`: script for systematically degrading images in a given dataset according to the degraded NA value `NA_L` based on the system parameters returned by `dataset_parameters`. Takes the following arguments:
    * `dataset` - (required) string, one of 'BreaKHis4X', 'BreaKHis10X', 'BreaKHis20X', 'BreaKHis40X', 'BACH', 'CAMELYON16-UMCU', 'CAMELYON16-RUMC', 'PCam'
    * `NA_L` - (required) positive float, should also be less than NA_H for this script to make sense but this is not currently enforced
    * `save_images` - (optional, default `true`) boolean flag of whether to save degraded images to disk
    * `sanity_check` - (optional, default `false`) boolean flag of whether to show PSF/MTF of system to ensure things are sensible
    * `plot_fft` - (optional, default `false`) boolean flag of whether to show image spectra of original/degraded images to ensure things are sensible
    * `Q` - (optional, default 10) integer, quantization/padding factor of pupil field size to ensure FFTs are computed correctly


## Installation
Tested on Windows 10 using MATLAB R2021a. Requires Image Processing Toolbox to be installed.

To cite this code, please use the following:
[![DOI](https://zenodo.org/badge/480455695.svg)](https://zenodo.org/badge/latestdoi/480455695)


