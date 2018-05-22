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
-D CMAKE_BUILD_TYPE=Debug ..

make
make test
make install
```

### current status:

see `examples/tastest.py`:

```
➜  examples git:(master) ✗ python tastest.py
TasmanianSG version: 6.0
TasmanianSG license: BSD 3-Clause with UT-Battelle disclaimer

          Grid Type:  Local Polynomial
         Dimensions:   2
            Outputs:   1
       Loaded nodes:   0
       Needed nodes:   145
               Rule:  Local polynomials
             Domain:  Canonical
              Order:   1
       Acceleration:  cpu-blas

```

julia:  

```
julia> Tasmanian.run()
6.0.0
Ptr{Void} @0x00007fcde8765270

          Grid Type:  Local Polynomial
         Dimensions:   2
            Outputs:   1
       Loaded nodes:   0
       Needed nodes:   145
               Rule:  Local polynomials
             Domain:  Canonical
              Order:   1
       Acceleration:  cpu-blas



```

