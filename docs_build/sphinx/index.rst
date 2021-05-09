.. rst-class:: emphasize-children

.. rst-class:: html-toggle

G1 and G2 fitting with clothoids, spline of clothods, circle arc and biarc
==========================================================================

**Github repository** https://github.com/ebertolazzi/Clothoids

**Matlab TOOLBOX** https://github.com/ebertolazzi/Clothoids/releases

Installation
------------


Matlab
~~~~~~

Download the toolbox and install as usual.
Run ``CompileClothoidsLib`` on Matlab console windows
to compile the ``mex`` of the library.

C++
~~~

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

Contents
--------

.. toctree::
    :maxdepth: 2

    Matlab_manual.rst
    api-cpp/root.rst
    api-matlab/root.rst


.. include:: authors.rst

.. include:: references.rst

License
~~~~~~~

.. literalinclude:: ../../License.txt
    :language: none
