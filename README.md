Clothoids [![Build Status](https://travis-ci.org/ebertolazzi/Clothoids.svg?branch=master)](https://travis-ci.org/ebertolazzi/Clothoids)
[![View ebertolazzi/Clothoids on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://it.mathworks.com/matlabcentral/fileexchange/64849-ebertolazzi-clothoids)

### G1 and G2 fitting with clothoids, spline of clothods, circle arc and biarc

**by Enrico Bertolazzi and Marco Frego**

for the documentation see `manual.md` or
[Doxygen documentation: http://ebertolazzi.github.io/Clothoids/](http://ebertolazzi.github.io/Clothoids/)

Installation
------------

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


Installation mex interface
--------------------------

The easy way is to download the compiled Toolbox at
[https://github.com/ebertolazzi/Clothoids/releases](https://github.com/ebertolazzi/Clothoids/releases) and install it.
After installation in Matlab run the command `CompileClothoidsLib`.

If you want to compile the Toolbox by yourself

~~~
cd toolbox
ruby populate_toolbox.rb
~~~

run Matlab and from the cmmand windows of MATLAB:

~~~
> cd toolbox
> setup
> CompileClothoidsLib
> open('../Clothoids.prj')
~~~

then compile toolbox and install.


**Authors:**
	
	Enrico Bertolazzi and Marco Frego
	Department of Industrial Engineering
	University of Trento
	enrico.bertolazzi@unitn.it
	m.fregox@gmail.com


References
----------

1. E. Bertolazzi, M. Frego,
   **G1 fitting with clothoids**,<br>
   Mathematical Methods in the Applied Sciences,
   John Wiley & Sons, 2014, vol. 38, n.5, pp. 881-897,<br>
   https://doi.org/10.1002/mma.3114

2. E. Bertolazzi, M. Frego,
   **On the G2 Hermite interpolation problem with clothoids**,<br>
   Journal of Computational and Applied Mathematics, 2018, vol. 15, n.341, pp. 99-116.<br>
   https://doi.org/10.1016/j.cam.2018.03.029

3. E. Bertolazzi, M. Frego,
   **Interpolating clothoid splines with curvature continuity**,<br>
   Mathematical Methods in the Applied Sciences, 2018, vol. 41, n.4, pp. 1099-1476.<br>
   https://doi.org/10.1002/mma.4700
   
4. E. Bertolazzi, M. Frego
   **A Note on Robust Biarc Computation**,<br>
   Computer-Aided Design & Applications 16 (5), 822-835
 

