## Medusa
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0).

Medusa is a C++14 compliant application, highly based on [HYDRA v.3](https://github.com/MultithreadCorner/Hydra), to perform physics data analysis of generic 4-body decays in massively parallel platforms on Linux systems. 


## Dependencies
Medusa depends on [HYDRA >= v.3.2.1](https://github.com/MultithreadCorner/Hydra), [GCC >= v.8](https://gcc.gnu.org/), [ROOT >= v.6.14](https://github.com/root-project/root), [libconfig >= v1.5](https://hyperrealm.github.io/libconfig/), [TCLAP >= v1.2.1](http://tclap.sourceforge.net/). Optionally  [CUDA >= 10.0](https://developer.nvidia.com/cuda-toolkit) is needed to use nVidia GPUs. 

## Installation, Build and Run the first tests
The first step is checkout [HYDRA v.3](https://github.com/MultithreadCorner/Hydra): 
```bash
mkdir <MedusaDevDir>
cd <MedusaDevDir>
git clone https://github.com/MultithreadCorner/Hydra.git Hydra
```
Then you can clone the Medusa repository:
```bash
git clone https://github.com/dbrundu/Medusa.git Medusa
```

Then you can setup the corresponding enveironment variables:
```bash
export CC=/usr/bin/gcc-8
export CXX=/usr/bin/g++-8
export HYDRA_INCLUDE_DIR=<path-to-hydra>
...
```

Starting from the Medusa folder, please create a `build` directory for convenience and run the cmake command:
```bash
cd Medusa
mkdir build
cd build
cmake -DHYDRA_INCLUDE_DIR=$HYDRA_INCLUDE_DIR ../

```

At this point you can run the first test (a preliminary benchmark for phis model): 
```bash
make test_phis_JpsiKK_tbb
./test_phis_JpsiKK_tbb -n=500000 -b=1
```
where `n` is the number of unewighted events to generate and `b` is a numerical boolean to perform the benchmark (1) or not (0).
If you have also OMP you could try to compile ```test_phis_JpsiKK_omp```.


