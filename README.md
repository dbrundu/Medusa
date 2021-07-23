## Medusa
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0).

Medusa is a C++14 compliant application, highly based on [HYDRA v.3](https://github.com/MultithreadCorner/Hydra), to perform physics data analysis of generic 4-body decays in massively parallel platforms on Linux systems. 


## Dependencies
Medusa depends on [HYDRA >= v.3.2.1](https://github.com/MultithreadCorner/Hydra), [GCC >= v.8](https://gcc.gnu.org/), [ROOT >= v.6.22.08](https://github.com/root-project/root), [libconfig >= v1.5](https://hyperrealm.github.io/libconfig/), [TCLAP >= v1.2.4](http://tclap.sourceforge.net/). Optionally  [CUDA >= v.9.2](https://developer.nvidia.com/cuda-toolkit) is needed to use nVidia GPUs. 
Medusa also uses [Catch2](https://github.com/catchorg/Catch2/tree/v2.x) to perform test and micro-benchmarks; this library is already included as unique header.

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

At this point you can run the first test (a preliminary benchmark for phis model with 1M events): 
```bash
make benchmarks_B0s_Jpsi_phi_tbb
./benchmarks_B0s_Jpsi_phi_tbb
```

You can run also the fit:
```bash
make fit_B0s_Jpsi_phi_tbb
./fit_B0s_Jpsi_phi_tbb
```



