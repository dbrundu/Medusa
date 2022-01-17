### CMakeList for Medusa examples and documentation
cmake_minimum_required(VERSION 3.2)

# project name
project(Medusa)

# warn user if system is not UNIX
if(NOT UNIX)
  message(FATAL_ERROR "This is an unsupported system.")
endif()

# cmake path dir
set(CMAKE_PREFIX_PATH "/home/aricci/Hydra")    # paths on cms server
set(MEDUSA_CMAKE_DIR "${PROJECT_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH "${MEDUSA_CMAKE_DIR}" ${CMAKE_MODULE_PATH})
set(CMAKE_VERBOSE_MAKEFILE  ON)

# check if compiler is C++14 compliant
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("--std=c++14" COMPILER_SUPPORTS_CXX14)
if(NOT COMPILER_SUPPORTS_CXX14)
  message(FATAL "The compiler ${CMAKE_CXX_COMPILER} has no C++14 support. Please use a different C++ compiler.")
endif()

# compiler flags
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
 MESSAGE(STATUS "Setting Clang flags")
 set(CMAKE_CXX_FLAGS " --std=c++14 -W -march=native -fPIC -O3 -ldl" CACHE STRING "compile flags" FORCE) 
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
 MESSAGE(STATUS "Setting GCC flags")
 set(CMAKE_CXX_FLAGS " --std=c++14 -W -march=native -fPIC -O4 -ldl" CACHE STRING "compile flags" FORCE)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
 MESSAGE(STATUS "Setting ICC flags")
 set(CMAKE_CXX_FLAGS " --std=c++14 -W -march=native -fPIC -O4 -ldl" CACHE STRING "compile flags" FORCE)
endif()

#-----------------------
# get Hydra
find_package(Hydra REQUIRED)
include_directories(${HYDRA_INCLUDE_DIR})
#-----------------------

#-----------------------
# get GSL
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
link_directories(${GSL_LIBRARIES})
#-----------------------

#-----------------------
# get TCLAP
find_package(TCLAP REQUIRED)
include_directories(${TCLAP_INCLUDE_DIR})
#-----------------------

#-----------------------
# get ROOT
find_package(ROOT COMPONENTS Minuit2)
if(ROOT_FOUND)
  include_directories(${ROOT_INCLUDE_DIR})
  link_directories(${ROOT_LIBRARY_DIR})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_ROOT_AVAILABLE_")
  if(${ROOT_Minuit2_LIBRARY} MATCHES "libMinuit2.so")
    set(Minuit2_FOUND ON)
  else(Minuit2_FOUND OFF)
  endif()
endif(ROOT_FOUND)
#-----------------------

#-----------------------
# get CUDA
find_package(CUDA 9.2)
if(CUDA_FOUND)
  link_directories( ${CUDA_TOOLKIT_ROOT_DIR}/lib64/)

  # set cuda flags
  SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}  -ftemplate-backtrace-limit=0; --cudart;  static; -O4; --expt-relaxed-constexpr; -fmad=true; --expt-extended-lambda;--relocatable-device-code=false ;
  -Xptxas -dlcm=ca; -Xptxas --opt-level=4)

  SET(CUDA_PROPAGATE_HOST_FLAGS ON)
  SET(CUDA_SEPARABLE_COMPILATION OFF)
  SET(CUDA_VERBOSE_BUILD OFF)

  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.8)
    LIST(APPEND CUDA_NVCC_FLAGS " -Xcompiler -D_MWAITXINTRIN_H_INCLUDED ")
  endif()

  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 5.4)
    LIST(APPEND CUDA_NVCC_FLAGS " -Xcompiler -D__CORRECT_ISO_CPP11_MATH_H_PROTO ")
  endif()

  # detect CUDA architecture and get best NVCC flags
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindCudaArch.cmake)
  SELECT_NVCC_ARCH_FLAGS(NVCC_FLAGS_EXTRA)
  LIST(APPEND CUDA_NVCC_FLAGS ${NVCC_FLAGS_EXTRA})
endif(CUDA_FOUND)
#-----------------------

#-----------------------
# get TBB
find_package(TBB)
if(TBB_FOUND)
  include_directories(${TBB_INCLUDE_DIRS})
  link_directories(${TBB_LIBRARY})
endif(TBB_FOUND)
#-----------------------

#-----------------------
# get OpenMP
find_package(OpenMP)
if(OPENMP_CXX_FOUND OR OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(OPENMP_CXX_FOUND OR OPENMP_FOUND)
#-----------------------

#-----------------------
# generate API documentation with Doxygen
find_package(Doxygen)
option(BUILD_DOCUMENTATION "Create and install the HTML based API documentation (requires Doxygen)" NO)
#-----------------------

#-----------------------
# get Medusa directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
#-----------------------

#-----------------------
# CPP builds always on
SET(BUILD_CPP_TARGETS  ON)

# crucial for Clang build, where openMP support is not embedded in the compiler
if( OPENMP_CXX_FOUND OR OPENMP_FOUND)
  SET(BUILD_OMP_TARGETS  ON)
else()
  SET(BUILD_OMP_TARGETS  OFF)
endif()

# not every computer has a NVIDIA GPU and/or CUDA installed 
if( CUDA_FOUND )
  SET(BUILD_CUDA_TARGETS  ON)
else()
  SET(BUILD_CUDA_TARGETS  OFF)
endif()

# enable TBB builds
if(TBB_FOUND)
  SET(BUILD_TBB_TARGETS  ON)
else()
  SET(BUILD_TBB_TARGETS  OFF)
endif(TBB_FOUND)
#-----------------------

#-----------------------
# messages
MESSAGE(STATUS "TCLAP Include path: ${TCLAP_INCLUDE_DIR}")
MESSAGE(STATUS "ROOT include path: ${ROOT_INCLUDE_DIR}")
MESSAGE(STATUS "ROOT library path: ${ROOT_LIBRARY_DIR}" )
MESSAGE(STATUS "ROOT libraries: ${ROOT_LIBRARIES}")
MESSAGE(STATUS "ROOT::Minuit2: ${ROOT_Minuit2_LIBRARY}")
MESSAGE(STATUS "Build CUDA/NVCC-based targets: ${BUILD_CUDA_TARGETS}")
MESSAGE(STATUS "CUDA include: ${CUDA_INCLUDE_DIRS}" )
MESSAGE(STATUS "CUDA libraries: ${CUDA_LIBRARIES}"  )
MESSAGE(STATUS "NVCC flags: ${CUDA_NVCC_FLAGS}" )
MESSAGE(STATUS "Build TBB targets: ${BUILD_TBB_TARGETS}")
MESSAGE(STATUS "TBB include: ${TBB_INCLUDE_DIRS}" )
MESSAGE(STATUS "TBB libraries: ${TBB_LIBRARIES}" )
MESSAGE(STATUS "Build OpenMP targets: ${BUILD_OMP_TARGETS}")
MESSAGE(STATUS "OpenMP flags: ${OpenMP_CXX_FLAGS}" )
MESSAGE(STATUS "OpenMP libraries:  ${OpenMP_CXX_LIBRARIES}" )
MESSAGE(STATUS "Host compiler flags: ${CMAKE_CXX_FLAGS}" )
MESSAGE(STATUS "Install directory: ${CMAKE_INSTALL_PREFIX}"  )
MESSAGE(STATUS "Project source dir: ${PROJECT_SOURCE_DIR}" )
MESSAGE(STATUS "Project build dir: ${PROJECT_BINARY_DIR}" )
#-----------------------

#+++++++++++++++++++++++++++
#   targets    +
#+++++++++++++++++++++++++++
include(${MEDUSA_CMAKE_DIR}/AddTarget_Fit.cmake)
include(${MEDUSA_CMAKE_DIR}/AddTarget_Benchmarks.cmake)

ADD_TARGET_FIT(fit_B0s_Jpsi_phi_signal)
ADD_TARGET_FIT(fit_B0s_Jpsi_phi_full)
ADD_TARGET_BENCHMARKS(benchmarks_complete_B0s_Jpsi_phi_signal)
ADD_TARGET_BENCHMARKS(benchmarks_simple_B0s_Jpsi_phi_signal)
ADD_TARGET_BENCHMARKS(benchmarks_complete_B0s_Jpsi_phi_full)
ADD_TARGET_BENCHMARKS(benchmarks_simple_B0s_Jpsi_phi_full)