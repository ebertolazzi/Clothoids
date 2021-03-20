Clothoids
=========

[![Build Status](https://travis-ci.org/ebertolazzi/Clothoids.svg?branch=master)](https://travis-ci.org/ebertolazzi/Clothoids)
[![View ebertolazzi/Clothoids on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://it.mathworks.com/matlabcentral/fileexchange/64849-ebertolazzi-clothoids)

G1 and G2 fitting with clothoids, spline of clothods, circle arc and biarc

*by Enrico Bertolazzi and Marco Frego*

for the documentation see
[online documentation](http://ebertolazzi.github.io/Clothoids/)

Installation
------------

Download the library

```sh
git clone git@github.com:ebertolazzi/Clothoids.git —recurse-submodules
```

if you forget ` —recurse-submodules` you must load the submdule quarticRootsFlocke next. To compile the library you can use `make`

```sh
make
```

or `cmake`

```sh
mkdir build
cd build
cmake ..
make
```

of `rake`

```sh
rake build_osx   # on mac
rake build_linux # on linux
rake build_win   # on windows
```


using `make` you have the library and headers in the following tree

```text
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
```

where OSTYPE can be `linux`, `osx`, `mingw_x64`, `win_x64`

### Online Documentation

Available at: [http://ebertolazzi.github.io/Clothoids](http://ebertolazzi.github.io/Clothoids)

**Authors:**
	
	Enrico Bertolazzi and Marco Frego
	Department of Industrial Engineering
	University of Trento
	enrico.bertolazzi@unitn.it
	m.fregox@gmail.com


References
----------

1. *E. Bertolazzi, M. Frego*,
   **G1 fitting with clothoids**,  
   Mathematical Methods in the Applied Sciences,  
   John Wiley & Sons, 2014, vol. 38, n.5, pp. 881-897,  
   [https://doi.org/10.1002/mma.3114]
   ([https://doi.org/10.1002/mma.3114)

2. *E. Bertolazzi, M. Frego*,
   **On the G2 Hermite interpolation problem with clothoids**,  
   Journal of Computational and Applied Mathematics,  
   2018, vol. 15, n.341, pp. 99-116.  
   [https://doi.org/10.1016/j.cam.2018.03.029]
   (https://doi.org/10.1016/j.cam.2018.03.029)

3. *E. Bertolazzi, M. Frego*,
   **Interpolating clothoid splines with curvature continuity**,  
   Mathematical Methods in the Applied Sciences,  
   2018, vol. 41, n.4, pp. 1099-1476.  
   [https://doi.org/10.1002/mma.4700](https://doi.org/10.1002/mma.4700)
   
4. *E. Bertolazzi, M. Frego*
   **A Note on Robust Biarc Computation**,  
   Computer-Aided Design & Applications 16 (5), 822-835  
   [http://www.cad-journal.net/files/vol_16/CAD_16(5)_2019_822-835.pdf]
   (http://www.cad-journal.net/files/vol_16/CAD_16(5)_2019_822-835.pdf)

