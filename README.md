# Monte-Carlo
Python script of Monte Carlo simulations for protein unfolding


- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Demo](#demo)
- [License](./LICENSE)
- [Citation](#citation)

# Overview
`MonteCarlo_DBM` is a Python script that conducts Monte Carlo simulation to study the dual-binding mode behavior for a single molecular system, that is capable of reproducing and validating the SMFS test (single molecule force sprectroscopy). Single moecule forced pulling test and corresponding force extension curves could be simulated based on the input parameters using a Bell-Evans model. 

The force could be loaded by either a constant pulling speed mode or a force clamp mode. Under constant speed mode, the cantilever is retracted with a constant pulling speed and the rupture force of the complex is recorded. Under force clamp mode, a fixed force is applied on the cantilever and the life time of the complex is recorded.


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

The versions of the software and packages are, specifically:
```
python        3.8.1
matplotlib    3.2.1
scipy         1.4.1
numpy         1.18.3
seaborn       0.10.1
```

# Installation Guide:

### Install from Github
```
git clone https://github.com/NashLab/Monte-Carlo
python3 setup.py install
```
- `sudo`, if required

Any missing packages in [Python Dependencies](#Python Dependencies) could be installed respectively using PIP installer:
```
pip install matplotlib
pip install scipy
pip install numpy
pip install seaborn
```

# Demo
An example output is given to simulate 1000 SMFS force curves at a pulling speed of 100 nm/s.
The python script is executed by the command:
```
python MonteCarlo_DBM.py
```

- Sampling rate: 12.5 points/nm
- Computation took 22.54 sec

All input parameters are included in `MonteCarlo_DBM.py`.

Results:
- 185/1000 curves showed low force (18.5%)
- 479/1000 curves showed high force without Xmod unfolding (47.9%)
- 336/1000 curves showed high force with Xmod unfolded (33.6%)
- Average Loading Rate of the complex rupture: 2456.546pN/s

Raw date could be found in [Demo](./Demo): 
- [Complex ruputure_force and Loading rate](./Demo/LR-RF_v100.0.txt)
- [Complex ruputure_histogram](./Demo/analysis-hist-ratio_0.2speed_100.0.pdf)
- [Xmod unfolding_force and Loading rate](./Demo/Xmod_v100.0.txt)
- [Xmod unfolding_histogram](./Demo/Xmod_unfolding-hist-ratio_0.2speed_100.0.pdf)
- [Logfile](./Demo/speed_100.0logfile.log)
 

# Citation
For usage of the package and associated manuscript, please cite according to the enclosed [citation.bib](./Demo/citation.bib).


