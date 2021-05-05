Clothoids
=========

|Build Status| |View ebertolazzi/Clothoids on File Exchange|

G1 and G2 fitting with clothoids, spline of clothoids, circle arc and
biarc

*by Enrico Bertolazzi (enrico.bertolazzi@unitn.it) and Marco Frego (marco.frego@unibz.it)*

for the documentation see `online documentation <http://ebertolazzi.github.io/Clothoids/>`__

Installation
------------

Download the library

.. code:: sh

   git clone git@github.com:ebertolazzi/Clothoids.git —recurse-submodules

if you forget ``—recurse-submodules`` you must load the submdule
quarticRootsFlocke next. To compile the library you can use ``make``

.. code:: sh

   make

or ``cmake``

.. code:: sh

   mkdir build
   cd build
   cmake ..
   make

of ``rake``

.. code:: sh

   rake build_osx   # on mac
   rake build_linux # on linux
   rake build_win   # on windows

using ``make`` you have the library and headers in the following tree

.. code:: text

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

where OSTYPE can be ``linux``, ``osx``, ``mingw_x64``, ``win_x64``

References
----------

1. | *E. Bertolazzi, M. Frego*, **G1 fitting with clothoids**,
   | Mathematical Methods in the Applied Sciences,
   | John Wiley & Sons, 2014, vol. 38, n.5, pp. 881-897,
   | `https://doi.org/10.1002/mma.3114 <https://doi.org/10.1002/mma.3114>`__

2. | *E. Bertolazzi, M. Frego*, **On the G2 Hermite interpolation
     problem with clothoids**,
   | Journal of Computational and Applied Mathematics,
   | 2018, vol. 15, n.341, pp. 99-116.
   | `https://doi.org/10.1016/j.cam.2018.03.029 <https://doi.org/10.1016/j.cam.2018.03.029>`__

3. | *E. Bertolazzi, M. Frego*, **Interpolating clothoid splines with
     curvature continuity**,
   | Mathematical Methods in the Applied Sciences,
   | 2018, vol. 41, n.4, pp. 1099-1476.
   | `https://doi.org/10.1002/mma.4700 <https://doi.org/10.1002/mma.4700>`__

4. | *E. Bertolazzi, M. Frego* **Point-Clothoid Distance and Projection Computation**,
   | SIAM Journal on Scientific Computing,
   | 2019, 41(5), A3326–A3353.
   | `https://doi.org/10.1137/18M1200439 <https://doi.org/10.1137/18M1200439>`__

5. | *E. Bertolazzi, M. Frego* **A Note on Robust Biarc Computation**,
   | Computer-Aided Design & Applications 16 (5), 822-835
   | `http://www.cad-journal.net/files/vol_16/CAD_16(5)_2019_822-835.pdf <http://www.cad-journal.net/files/vol_16/CAD_16(5)_2019_822-835.pdf>`__

.. |Build Status| image:: https://travis-ci.org/ebertolazzi/Clothoids.svg?branch=master
   :target: https://travis-ci.org/ebertolazzi/Clothoids
.. |View ebertolazzi/Clothoids on File Exchange| image:: https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg
   :target: https://it.mathworks.com/matlabcentral/fileexchange/64849-ebertolazzi-clothoids
