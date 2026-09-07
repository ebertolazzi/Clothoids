# Clothoids: Getting Started {#mainpage}

**Clothoids** is a C++20 library for building, evaluating, and fitting curves
made of straight segments, circular arcs, biarcs, clothoids (Euler spirals),
piecewise-clothoid lists, and Dubins paths. It is the API reference generated
directly from the header comments; for build instructions and the Python/MATLAB
bindings see the project [README](https://github.com/ebertolazzi/Clothoids).

## The curve hierarchy

Every curve type in the library derives from G2lib::BaseCurve, which defines
the common interface (`length()`, `theta()`, `X()`, `Y()`, offset curves,
bounding boxes, intersections, ...). The concrete curve types are:

| Class                    | Curve type                                   |
|---------------------------|-----------------------------------------------|
| G2lib::LineSegment        | straight line segment                          |
| G2lib::CircleArc          | circular arc                                   |
| G2lib::Biarc              | a pair of tangent circular arcs (G1 biarc)     |
| G2lib::BiarcList          | a list of consecutive biarcs                   |
| G2lib::ClothoidCurve       | a single clothoid (Euler spiral) arc           |
| G2lib::ClothoidList       | a list of consecutive clothoid arcs (G2 spline)|
| G2lib::PolyLine           | a polyline (list of line segments)             |
| G2lib::Dubins              | a 2-point Dubins path (3 arcs)                 |
| G2lib::Dubins3p            | a 3-point Dubins path                          |

G2lib::CurveType (`G2lib::CurveType::LINE`, `CIRCLE`, `BIARC`, `BIARC_LIST`,
`CLOTHOID`, `CLOTHOID_LIST`, `DUBINS`, `DUBINS3P`, `POLYLINE`) identifies the
runtime type of a curve reached through a `G2lib::BaseCurve` pointer or
reference; see G2lib::BaseCurve::type().

## Minimal example: building a clothoid

A clothoid arc can be built directly from its standard parameters
(initial position/angle, curvature, curvature rate, length):

```cpp
#include "Clothoids.hh"

using G2lib::real_type;
using Utils::m_pi;

int main() {
  G2lib::ClothoidCurve clothoid( "my_clothoid" );
  // x0, y0, theta0, kappa0, dk, length
  clothoid.build( 0.0, 0.0, m_pi / 4, 0.0, 0.1, 10.0 );

  real_type const L  = clothoid.length();
  real_type const xm = clothoid.X( L / 2 );
  real_type const ym = clothoid.Y( L / 2 );
  real_type const th = clothoid.theta( L / 2 );
}
```

## Fitting a clothoid between two oriented points (\f$G^1\f$ Hermite problem)

The library's signature feature is solving for the clothoid (or list of
clothoids) that interpolates two points with prescribed tangent directions.
See G2lib::ClothoidCurve::build_G1() for the single-arc case:

```cpp
G2lib::ClothoidCurve clothoid( "G1" );
int iter = clothoid.build_G1(
  0.0, 0.0, m_pi / 6,   // x0, y0, theta0
  10.0, 4.0, -m_pi / 6  // x1, y1, theta1
);
```

For a smooth multi-arc spline through an ordered sequence of waypoints, see
G2lib::ClothoidList and its `build*` / `G2solve*` family of methods.

## Where to go next

- G2lib::BaseCurve — the common curve interface implemented by every class above.
- G2lib::ClothoidCurve, G2lib::ClothoidList — clothoid arcs and clothoid splines.
- G2lib::Dubins, G2lib::Dubins3p — shortest paths with a bounded turning radius.
- G2lib::Biarc, G2lib::BiarcList — circular-arc alternative to clothoid splines.
- The [namespace G2lib](namespaces_namespaces.html) index lists every free
  function (collision detection, intersections, curve-type promotion, ...).

*This page is generated from `doc_generate/tutorial.md` (Doxygen mainpage).
See also the project's Sphinx-based documentation for narrative guides.*
