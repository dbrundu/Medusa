## Medusa
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0).

Medusa is a C++14 compliant application, highly based on [HYDRA v.3](https://github.com/MultithreadCorner/Hydra), designed to perform physics data analyses of generic 4-body decays deploying massively parallel platforms (TBB, OpenMP, CUDA) on Linux systems.


## Dependencies
Medusa depends on [HYDRA >= v.3.2.1](https://github.com/MultithreadCorner/Hydra), [GCC >= v.8](https://gcc.gnu.org/), [ROOT >= v.6.22.08](https://github.com/root-project/root), [libconfig >= v1.5](https://hyperrealm.github.io/libconfig/), [TCLAP >= v1.2.4](http://tclap.sourceforge.net/). Optionally  [CUDA >= v.9.2](https://developer.nvidia.com/cuda-toolkit) is needed to use nVidia GPUs. 
Medusa also uses [Catch2](https://github.com/catchorg/Catch2/tree/v2.x) to perform test and micro-benchmarks; this library is already included as unique header.

## Installation, Build and Run the first tests
1- Install all dependences, as ROOT and CMAKE.

2- Clone [HYDRA v.3](https://github.com/MultithreadCorner/Hydra):
```bash
mkdir <MedusaDevDir>
cd <MedusaDevDir>
git clone https://github.com/MultithreadCorner/Hydra.git Hydra
```
3- Clone the Medusa repository:
```bash
git clone https://github.com/dbrundu/Medusa.git Medusa
```

4- Setup the corresponding environment variables:
```bash
export CC=/usr/bin/gcc-8
export CXX=/usr/bin/g++-8
export HYDRA_INCLUDE_DIR=<path-to-hydra>
...
```

5- Starting from the Medusa folder, please create a `build` directory for convenience and run the commands:
```bash
cd Medusa
mkdir build
cd build
cmake -DHYDRA_INCLUDE_DIR=$HYDRA_INCLUDE_DIR ../
```

5.1- Alteratively, open with an editor the file /Medusa/CMakeLists.txt and under the line "# cmake path dir" insert the command
```bash
set(CMAKE_PREFIX_PATH "<path-to-hydra>")
```
Then,
```bash
cd Medusa
mkdir build
cd build
cmake ../
```

6- To compile all, run the command
```bash
make
```

7- Run the executables, for example
```bash
./fit_sim_phis_full_model_tbb
```



