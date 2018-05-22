# Tasmanian [WIP]

[![Build Status](https://travis-ci.com/floswald/Tasmanian.jl.svg?branch=master)](https://travis-ci.com/floswald/Tasmanian.jl)

wrapping [tasmanian](https://github.com/ORNL/Tasmanian)

## Installation

To install julia package

```julia
Pkg.clone("git@github.com:floswald/Tasmanian.jl.git")
```

in terminal, to install the Tasmanian library.
this will install the library into dir `/Applications/TSG`


```
git clone https://github.com/ORNL/TASMANIAN
cd TASMANIAN
mkdir build
cd build
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
make
make test
make install
```

### use as python module

see `examples/tastest.py`:

```
➜  examples git:(master) ✗ python tastest.py 
TasmanianSG version: 6.0
TasmanianSG license: BSD 3-Clause with UT-Battelle disclaimer

-------------------------------------------------------------------------------------------------
Example 1 for OSM: interpolate f(x,y) = cos(0.5 * pi * x) * cos(0.5 * pi * y)
       using fixed sparse grid with depth 5
       the error is estimated as the maximum from 1000 random points

 For localp    Number of points: 145   Max. Error: 7.9623370348109734e-03
```

