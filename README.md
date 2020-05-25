## Medusa
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

If you want to run a first test you need to checkout the Hydra library v.3, currently under development. Please create a Hydra folder somewhere and clone the repository:
```bash
git clone --single-branch --branch develop-hydra3 https://github.com/MultithreadCorner/Hydra.git
```
Then you can setup the corresponding enveironment variable and compile the Medusa test application in the following way (starting from the Medusa folder, create a build directory for convenience):
```bash
mkdir build
cd build
export HYDRA_INCLUDE_DIR=<path-to-hydra>
cmake -DHYDRA_INCLUDE_DIR=$HYDRA_INCLUDE_DIR ../
make test_cpp
./test_cpp
```
If you have TBB or OMP try also to compile ```test_tbb``` or ```test_omp```.
