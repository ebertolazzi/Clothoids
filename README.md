Clothoids
=========

[![](https://travis-ci.org/ebertolazzi/Clothoids.svg?branch=master)](https://travis-ci.org/ebertolazzi/Splines) [![](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://it.mathworks.com/matlabcentral/fileexchange/64849-ebertolazzi-clothoids)

G1 and G2 fitting with clothoids, spline of clothoids, circle arc and
biarc.

This library contains the implementation of the algorithms on clothoids, splines of clothoids,  arc, bi-arc splines of biarcs described in the works:

- **E.Bertolazzi, M.Frego**,
  G1 fitting with clothoids,
  Mathematical Methods in the Applied Sciences, 2015,
  https://doi.org/10.1002/mma.3114

- **E.Bertolazzi, M.Frego**,
  Interpolating clothoid splines with curvature continuity,
  Mathematical Methods in the Applied Sciences, 2018,
  https://doi.org/10.1002/mma.4700

- **E.Bertolazzi, M.Frego**,
  On the G2 Hermite interpolation problem with clothoids,
  Journal of Computational and Applied Mathematics, 2018,
  https://doi.org/10.1016/j.cam.2018.03.029

- **E.Bertolazzi, M.Frego**,
  A Note on Robust Biarc Computation,
  Computer-Aided Design And Applications, 2019,
  https://doi.org/10.14733/cadaps.2019.822-835

- **E.Bertolazzi, M.Frego**,
  Point-Clothoid Distance and Projection Computation,
  SIAM Journal on Scientific Computing, 2019,
  https://doi.org/10.1137/18M1200439

- **E.Bertolazzi, M.Frego, Francesco Biral**,
  Interpolating splines of biarcs from a sequence of planar points,
  Computer-Aided Design And Applications, 2020,
  https://doi.org/10.14733/cadaps.2021.66-85

- **E.Bertolazzi, P.Bevilacqua, M.Frego**,
  Efficient intersection between splines of clothoids,
  Mathematics and Computers in Simulation, 2020,
  https://doi.org/10.1016/j.matcom.2019.10.001

- **E.Bertolazzi, M.Frego, Francesco Biral**,
  Interpolating splines of biarcs from a sequence of planar points,
  Computer-Aided Design And Applications, 2020,
  https://doi.org/10.14733/cadaps.2021.66-85

- **M.Frego**,
  Closed form parametrisation of 3D clothoids by arclength
with both linear varying curvature and torsion,
  Applied Mathematics and Computation, 2022,
  https://doi.org/10.1016/j.amc.2021.126907

A clothoid is a curve described by the parametric eqautions:

$$
  x(s)=\int_0^s \cos\left(\frac{1}{2}\kappa'\tau^2+\kappa_0\tau+\vartheta_0\right)\,\mathrm{d}\tau
$$

$$
  y(s)=\int_0^s \sin\left(\frac{1}{2}\kappa'\tau^2+\kappa_0\tau+\vartheta_0\right)\,\mathrm{d}\tau
$$

when $\kappa'=0$ the clothoids reduce to a circle arc and when $\kappa'=\kappa_0=0$ the clothoids reduce to a straight segment.

The library contains the following objects:

- Segment
- Circle Arc
- Clothoids
- Bi-arc
- spline of
  - Segment
  - Circle Arc
  - Clothoids (with G1 and G2 continuity)
  - Bi-arc
- triangles
- bounding box

and fast algorithms involving the objects, in particular:

- evaluation
- intersection (between objects)
- point-object distance

Library is written in `C++11` with a `MATLAB` mex interface. Thus can be used in fast compiled application or in `MATLAB` scripts.

To compile the `C++11` library the easy way require `cmake` and `rake`

```
ruby setup.rb
```

then

```
rake
```

for more details see: **online documentation** at http://ebertolazzi.github.io/Clothoids/
