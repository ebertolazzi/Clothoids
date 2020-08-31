# Python Wrapper

The python wrapper is one to one with the Clothoids library. 

## Requirements

The wrapper uses the Github version of `pybind11`, and it is not compatible with 
the `pip` version (that lacks cmake support).

In order to compile `pybind11` (commands valid on a Unix machine, Win version should
apply some minor modifications):

``` bash
export PYTHON_TARGET=/usr/bin/python3 # or /usr/bin/python2 for py2 version
git clone https://github.com/pybind/pybind11.git
cd pybind11
mkdir build && cd build
cmake -DPYTHON_EXECUTABLE:FILEPATH=$PYTHON_TARGET ..
sudo make install -j8
```

## Compile the wrapper

It is possible to compile the wrapper for both python 2 and 3. In order to compile
the wrapper:

``` bash
export PYTHON_TARGET=/usr/bin/python3 # or /usr/bin/python2 for py2 version
cd src_py
mkdir build && cd build
cmake -DPYTHON_EXECUTABLE:FILEPATH=$PYTHON_TARGET ..
make install -j8
```

The shared object to be imported in python is stored in `{clothoid_dir}/lib/lib/G2lib.{ext}`,
near the Clothoids static library.

## Usage

In python:

``` python
import sys
sys.path.insert(0, "{clothoid_dir}/lib/lib")

import G2lib

help(G2lib)  # This is the best way to know the actual python interface
```

The wrapper wraps one to one the original C++ library, almost the same interface should be expected.

## Using Levemberg-Marquardt (Eigen3) as interpolator

If Eigen3 is installed via `apt`, the library will append to `G2lib` a function `interpolator` that
solves the interpolation optimal problem. Required pckages are:

``` bash
sudo apt install libeigen3-dev
```

This enable `G2lib.buildP[1-2]`


## Using IPOPT as interpolator

If IPOPT is installed via `apt`, the library will append to `G2lib` a function `interpolator` that
solves the interpolation optimal problem. Required pckages are:

``` bash
sudo apt install coinor-libipopt-dev coinor-libipopt1v5
```

This enable `G1lib.buildP[3-9]`

## TODO

 - [x] Include in Clothoids repo
 - [x] Implement an interpolator in Ipopt
 - [ ] Add docstrings to all methods
 - [ ] Naming arguments in methods