# Monte-Carlo
Monte Carlo simulations of protein unfolding


- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Setting up the development environment](#setting-up-the-development-environment)
- [License](./LICENSE)
- [Citation](#citation)

# Overview
`MonteCarlo_DBM` is a Python script that conducts Monte Carlo simulation of the dual-binding mode behavior for a single molecular system, that is capable of reproducing and validating the SMFS test (single molecule force sprectroscopy).


# System Requirements
## Hardware requirements
`MonteCarlo_DBM` python script requires only a standard computer with enough RAM to support the in-memory operations.

## Software requirements
### OS Requirements
This package is supported for *Windows*, *macOS* and *Linux*. The package has been tested on the following systems:
+ Windows:  Windows 10.1909
+ macOS:    Mojave  10.14.1
+ Linux:    Ubuntu  16.04

### Python Dependencies
`MonteCarlo_DBM` mainly depends on the Python scientific stack.
```
matplotlib
scipy
numpy
seaborn
```

# Installation Guide:

### Install from Github
```
git clone https://github.com/NashLab/Monte-Carlo
python3 setup.py install
```
- `sudo`, if required

# Setting up the development environment:

# Citation
For usage of the package and associated manuscript, please cite according to the enclosed [citation.bib](./citation.bib).
