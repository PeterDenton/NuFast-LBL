| arXiv |
|:-----:|
|[![arXiv](https://img.shields.io/badge/arXiv-2405.02400-orange.svg)](https://arXiv.org/abs/2405.02400)|

# **NuFast-LBL**: A fast code for long-baseline neutrino oscillation probabilities in matter

### Overview
**NuFast-LBL** is designed to calculate all nine neutrino oscillation probabilities in matter for long-baseline accelerator (e.g. NOvA, T2K, DUNE, HK) and reactor experiments (e.g. JUNO) very quickly, using the algorithms optimized for realistic oscillation scenarios. **NuFast-LBL** is provided in fortran, c++, and python, although no particular guarantees are made that the python code is "fast".

### General usage
The folders for each language contain a minimal working implementation to calculate the probabilities in matter and vacuum. The functions take as input the six oscillation parameters, with the mixing angles in the form of `s12sq` and so on. The CP violating phase, `delta`, is in radians. The atmospheric mass ordering is set by the sign of `Dmsq31`. The mass-squared differences are in eV<sup>2</sup>. For antineutrinos change only the energy to negative. The neutrino energy is in GeV, the baseline is in km, and the density is in g/cc. The electron fraction, `Ye`, is typically about 0.5 in the Earth. The function returns a 3x3 array `probs_returned` which contains all nine channels. The first index corresponds to the initial neutrino flavor (fortran is 1-indexed so e=1, while c++ and python are 0-indexed so e=0). Finally, the `N_Newton` parameter provides a proxy for the precision. `N_Newton=0` is already precise enough for even the next generation of neutrino oscillation experiments and higher values add orders of magnitude to the precision for minimal computational cost. Going beyond 2 is never necessary. If you are using unusual oscillation parameters (e.g. `|Dmsq31|<Dmsq21`), use a negative number which calls the (slower) exact analytic expression for `lambda3`.

We strongly encourage a careful examination of the constants used (G<sub>F</sub>, etc.) when comparing to existing codes.

An additional function, `Probability_Vacuum_LBL`, is provided for vacuum calculations. The algorithm is the same as the standard one, but it skips the matter effect.

The code can be compiled with the simple `compile.sh` file provided in the fortran and c++ folders. The benchmark codes can be compiled with the given makefiles by running `make`. The inclusion of aggressive compiler flags like `-ffast-math` can be included or not.

If a user requires a python code, we recommend the use of cython, f2py, or something similar from the c++ or fortran codes. The python code is included so that those who are only familiar with python can easily parse the algorithm itself, not for efficient computations.

### Benchmark
The benchmark folder contains the same code in a slightly different structure, as well as some other neutrino oscillation probability codes. It computes the precision of the **NuFast-LBL** algorithm and also benchmarks the computational speed.

### GLoBES
The included `GLoBES` code is adapted from the well known GLoBES code at [https://www.mpi-hd.mpg.de/personalhomes/globes/index.html](https://www.mpi-hd.mpg.de/personalhomes/globes/index.html). If that portion of the code is used, please cite them as appropriate.

Running the GLoBES code requires GSL. The makefile should work fine if it is installed in a typical Linux system.

### Usage
If you use this code, please cite the associated paper [arXiv:2405.02400](https://arxiv.org/abs/2405.02400) by Peter Denton and Stephen Parke. Please also let us know if you find any bugs or further optimizations or if you run your own speed tests.
