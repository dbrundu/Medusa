## Medusa
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

If you want to run a first test you need to checkout the Hydra library v.3, recently released and available on the master branch. Please create a Hydra folder somewhere and clone the repository:
```bash
git clone https://github.com/MultithreadCorner/Hydra.git
```
Then you can setup the corresponding enveironment variable and compile the Medusa test application (a preliminary benchmark for phis model) in the following way: starting from the Medusa folder, please create a build directory for convenience:
```bash
mkdir build
cd build
export HYDRA_INCLUDE_DIR=<path-to-hydra>
cmake -DHYDRA_INCLUDE_DIR=$HYDRA_INCLUDE_DIR ../
make test_phis_JpsiKK_tbb
./test_phis_JpsiKK_tbb -n=500000 -b=1
```
where n is the number of unewighted events to generate and b is a numerical boolean to perform the benchmark (1) or not (0).
If you have also OMP you could try to compile ```test_phis_JpsiKK_omp```.
