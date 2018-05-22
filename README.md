# Tasmanian

[![Build Status](https://travis-ci.com/floswald/Tasmanian.jl.svg?branch=master)](https://travis-ci.com/floswald/Tasmanian.jl)

wrapping [tasmanian](https://github.com/ORNL/Tasmanian)

## how to

* build the library for a local system? (i.e. travis?)

```bash
cmake \
-D CMAKE_INSTALL_PREFIX=/Applications/TSG \
-D Tasmanian_STRICT_OPTIONS=OFF \
-D Tasmanian_ENABLE_BLAS=ON \
-D Tasmanian_ENABLE_PYTHON=ON \
-D Tasmanian_ENABLE_MATLAB=OFF \
-D Tasmanian_ENABLE_FORTRAN=ON \
-D Tasmanian_ENABLE_CUBLAS=OFF \
-D Tasmanian_ENABLE_CUDA=OFF \
-D Tasmanian_ENABLE_MPI=OFF \
-D Tasmanian_SHARED_LIBRARY=ON \
-D Tasmanian_STATIC_LIBRARY=ON \
-D CMAKE_BUILD_TYPE=Debug \
..
```