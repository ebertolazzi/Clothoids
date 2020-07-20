Clothoids [![Build Status](https://travis-ci.org/ebertolazzi/Clothoids.svg?branch=master)](https://travis-ci.org/ebertolazzi/Clothoids)

### G1 and G2 fitting with clothoids, spline of clothods, circle arc and biarc

**by Enrico Bertolazzi and Marco Frego**

for the documentation see `manual.md` or
[Doxygen documentation: http://ebertolazzi.github.io/Clothoids/](http://ebertolazzi.github.io/Clothoids/)

**installation**

download the library

~~~~
git clone git@github.com:ebertolazzi/Clothoids.git —recurse-submodules
~~~~

if you forget ` —recurse-submodules` you must load the submdule quarticRootsFlocke next.

to compile the library you can use `make`

~~~~
make
~~~~

or `cmake`

~~~~
mkdir build
cd build
cmake ..
make
~~~~

using `make` you have the library and headers in the following tree

~~~~
`-- lib
    |-- include
    |   |-- AABBtree.hh
    |   |-- Biarc.hh
    |   |-- BiarcList.hh
    |   |-- Circle.hh
    |   |-- Clothoid.hh
    |   |-- ClothoidAsyPlot.hh
    |   |-- ClothoidList.hh
    |   |-- Fresnel.hh
    |   |-- G2lib.hh
    |   |-- Line.hh
    |   |-- PolyLine.hh
    |   |-- PolynomialRoots-Utils.hh
    |   |-- PolynomialRoots.hh
    |   `-- Triangle2D.hh
    `-- lib
        |-- libClothoids_OSTYPE.dylib
        `-- libClothoids_OSTYPE_static.a
~~~~

where OSTYPE can be `linux`, `osx`, `mingw_x64`, `win_x64`


**Authors:**
	
	Enrico Bertolazzi and Marco Frego
	Department of Industrial Engineering
	University of Trento
	enrico.bertolazzi@unitn.it
	m.fregox@gmail.com

