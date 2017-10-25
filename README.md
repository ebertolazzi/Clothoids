G1fitting [![Build Status](https://travis-ci.org/ebertolazzi/G1fitting.svg?branch=master)](https://travis-ci.org/ebertolazzi/G1fitting)

### G1 and G2 fitting with clothoids
**by Enrico Bertolazzi and Marco Frego**

The script `buildClothoid` implements the algorithm described in the paper

*G1 fitting with clothoids*, Mathematical Methods in the Applied Sciences, John Wiley & Sons, (2014), Ltd,.
http://onlinelibrary.wiley.com/doi/10.1002/mma.3114/abstract

**Description:**
Given two points and two direction associated with the points, 
a clothoid, i.e. a curve with linear varying curvature is computed
in such a way it pass to the points with the prescribed direction.
The solution in general is not unique but chosing the one for
which the angle direction variation is less than `2*pi` the solution
is unique.
The sofware solves the nonlinear system associated to the fitting problem
computing initial curvature and its derivative with the lenght of the curve.
An additional routine for the computation of the points along a clothoid
curve is added for convenience.

**Usage:**
To compute curvature and length of a clothoid passing throught points
P0=`(x0,y0)`, P1=`(x1,y1)` with angles `theta0` and `theta1` use the
function `buildClothoid`

`clot = buildClothoid( x0, y0, theta0, x1, y1, theta1, tol ) ;`

The parameter `tol` (usually `1e-10`) is a tolerance parameter
used to stop Newton iteration.
The resulting curve can be described by the 5 parameters

  - `(clot.x0, clot.y0)` initial point
  - `clot.theta0`        initial direction (angle)
  - `clot.k0`            initial curvature of the curve
  - `clot.dk`            derivative of the cuvature along arc length

plus a 6th parameter `clot.L`

  - `clot.L`         total length of the curve connecting P0 and P1

to compute points along a clothoid curve use the function `pointsOnClothoid`

`XY = pointsOnClothoid( x0, y0, theta0, k, dk, 0:L/npts:L ) ;

or

`XY = pointsOnClothoid( clot, 0: clot.L/npts:clot.L ) ;`

This function uses the 5 parameters `x0`, `y0`, `theta0`, `k`, `dk`
which indentify the curve. The parameter `L` is used to determine length 
of the portion of the curve to compute. The parameter `npts` is the number
of points computed along the curve. 

`XY` is a `2 x npts` matrix whose columns are the points along the curve. 

To plot the computed curve use MATLAB `plot` command as usual:

`plot( XY(1,:), XY(2,:), '-r' ) ;`

Four sample scripts: TestN0, TestN1, TestN2 and TestN3 shows how to use the functions.

**Mex files for fast computation:**

In directory `src_mex` you find a C++ implementation of the proposed algorithm 
with `mex` interface. To compile run `Compile` from MATLAB window.
After compilation the compiled version of the scripts

- buildClothoid
- evalClothoid
- FresnelCS
- GeneralizedFresnelCS
- pointsOnClothoid

are available in the `matlab` directory.

**Additional mex**

The mex implementation of the script:

- intersectClothoid
- TriTriOverlap
- biarc

are added to the library. The first compute all the intersections between two clothoids. The second check if two triangles (planar) overlap (used in the intersection computation).
The third compute a biarc given G1 data.

The experimental script

- buildClothoid2arcG2
- buildClothoid3arcG2

implements the algorithms described in the paper
*On the G2 Hermite Interpolation Problem with Clothoids*, submitted for publication.

**Authors:**
	
	Enrico Bertolazzi and Marco Frego
	Department of Industrial Engineering
	University of Trento
	enrico.bertolazzi@unitn.it
	m.fregox@gmail.com

