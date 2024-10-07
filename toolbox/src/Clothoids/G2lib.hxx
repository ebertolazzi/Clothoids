/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

///
/// file: G2lib.hh
///

//!
//! Clothoid computations routine namespace
//!
//! @html_image{Art_Figure_1.png,width=60%}
//!
namespace G2lib {

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  extern real_type const m_1_sqrt_pi;  //!< \f$ 1/\sqrt{\pi} \f$
  extern real_type const machepsi;     //!< machine espilon \f$ \varepsilon \f$
  extern real_type const machepsi10;   //!< \f$ 10\varepsilon \f$
  extern real_type const machepsi100;  //!< \f$ 100\varepsilon \f$
  extern real_type const machepsi1000; //!< \f$ 1000\varepsilon \f$
  extern real_type const sqrtMachepsi; //!< \f$ \sqrt{\varepsilon} \f$
  extern bool            intersect_with_AABBtree;

  extern integer const  G2LIB_AABB_CUT;
  extern integer const  G2LIB_AABB_MIN_NODES;

  // for CLOTHOIDS_BACK_COMPATIBILITY
  extern bool use_ISO;

  #endif

  #ifdef CLOTHOIDS_BACK_COMPATIBILITY

  static
  inline
  void
  lib_use_ISO()
  { use_ISO = true; }

  static
  inline
  void
  lib_use_SAE()
  { use_ISO = false; }

  #endif

  //!
  //! Disable AABB tree in computation
  //!
  static
  inline
  void
  noAABBtree()
  { intersect_with_AABBtree = false; }

  //!
  //! Enable AABB tree in computation
  //!
  static
  inline
  void
  yesAABBtree()
  { intersect_with_AABBtree = true; }

  /*
   * sin(x)/x
   */
  real_type Sinc( real_type x );     //!< \f$ \frac{\sin x}{x} \f$
  real_type Sinc_D( real_type x );   //!< \f$ \frac{\mathrm{d}}{\mathrm{d}x} \frac{\sin x}{x} \f$
  real_type Sinc_DD( real_type x );  //!< \f$ \left(\frac{\mathrm{d}}{\mathrm{d}x}\right)^2 \frac{\sin x}{x} \f$
  real_type Sinc_DDD( real_type x ); //!< \f$ \left(\frac{\mathrm{d}}{\mathrm{d}x}\right)^3 \frac{\sin x}{x} \f$

  /*
   * (1-cos(x))/x
   */
  real_type Cosc( real_type x );     //!< \f$ \frac{1-\cos x}{x} \f$
  real_type Cosc_D( real_type x );   //!< \f$ \frac{\mathrm{d}}{\mathrm{d}x} \frac{1-\cos x}{x} \f$
  real_type Cosc_DD( real_type x );  //!< \f$ \left(\frac{\mathrm{d}}{\mathrm{d}x}\right)^2 \frac{1-\cos x}{x} \f$
  real_type Cosc_DDD( real_type x ); //!< \f$ \left(\frac{\mathrm{d}}{\mathrm{d}x}\right)^3 \frac{1-\cos x}{x} \f$

  /*
   * atan(x)/x
   */
  real_type Atanc( real_type x );     //!< \f$ \frac{\arctan x}{x} \f$
  real_type Atanc_D( real_type x );   //!< \f$ \frac{\mathrm{d}}{\mathrm{d}x} \frac{\arctan x}{x} \f$
  real_type Atanc_DD( real_type x );  //!< \f$ \left(\frac{\mathrm{d}}{\mathrm{d}x}\right)^2 \frac{\arctan x}{x} \f$
  real_type Atanc_DDD( real_type x ); //!< \f$ \left(\frac{\mathrm{d}}{\mathrm{d}x}\right)^3 \frac{\arctan x}{x} \f$

  //!
  //! Add or remove multiple of \f$ 2\pi \f$ to an angle  in order to put it in the range \f$ [-\pi,\pi]\f$.
  //!
  void rangeSymm( real_type & ang );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //!
  //! Return minumum and maximum of three numbers
  //!
  inline
  void
  minmax3(
    real_type   a,
    real_type   b,
    real_type   c,
    real_type & vmin,
    real_type & vmax
  ) {
    vmin = vmax = a;
    if ( b < vmin ) vmin = b;
    else            vmax = b;
    if      ( c < vmin ) vmin = c;
    else if ( c > vmax ) vmax = c;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //!
  //! Project point `(qx,qy)` to the circle arc passing from `(x0,y0)`
  //! with tangent direction `(c0,s0)` curvature `k` length `L`
  //!
  //! \param[in] x0 \f$ x \f$-starting point of circle arc
  //! \param[in] y0 \f$ y \f$-starting point of circle arc
  //! \param[in] c0 \f$ \cos \theta_0 \f$
  //! \param[in] s0 \f$ \sin \theta_0 \f$
  //! \param[in] k  curvature
  //! \param[in] L  arc length
  //! \param[in] qx \f$ x \f$-point to be projected
  //! \param[in] qy \f$ y \f$-point to be projected
  //!
  //! \return distance point circle
  //!

  real_type
  projectPointOnCircleArc(
    real_type x0,
    real_type y0,
    real_type c0,
    real_type s0,
    real_type k,
    real_type L,
    real_type qx,
    real_type qy
  );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //!
  //! Project point `(qx,qy)` to the circle passing from `(x0,y0)`
  //! with tangent direction `(c0,s0)` and curvature `k`
  //!
  //! \param[in] x0     x-starting point of circle arc
  //! \param[in] y0     y-starting point of circle arc
  //! \param[in] theta0 initial angle
  //! \param[in] k      curvature
  //! \param[in] qx     x-point to be projected
  //! \param[in] qy     y-point to be projected
  //!
  //! \return distance point circle
  //!

  real_type
  projectPointOnCircle(
    real_type x0,
    real_type y0,
    real_type theta0,
    real_type k,
    real_type qx,
    real_type qy
  );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //!
  //! Check if point `(qx,qy)` is inside the circle passing from `(x0,y0)`
  //! with tangent direction `(c0,s0)` and curvature `k`
  //!
  //! \param[in] x0 starting \f$x\f$-coordinate of the circle arc
  //! \param[in] y0 starting \f$y\f$-coordinate of the circle arc
  //! \param[in] c0 \f$ \cos \theta \f$
  //! \param[in] s0 \f$ \sin \theta \f$
  //! \param[in] k  Curvature of the circle
  //! \param[in] qx \f$ x \f$-coordinate point to check
  //! \param[in] qy \f$ y \f$-coordinate point to check
  //!
  //! \return true if point is inside
  //!
  inline
  bool
  pointInsideCircle(
    real_type x0,
    real_type y0,
    real_type c0,
    real_type s0,
    real_type k,
    real_type qx,
    real_type qy
  ) {
    real_type cx  = x0 - s0/k;
    real_type cy  = y0 + c0/k;
    real_type dst = hypot( qx - cx, qy - cy );
    return dst*k <= 1;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  solveLinearQuadratic(
    real_type A,
    real_type B,
    real_type C,
    real_type a,
    real_type b,
    real_type c,
    real_type x[],
    real_type y[]
  );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  solveLinearQuadratic2(
    real_type A,
    real_type B,
    real_type C,
    real_type x[],
    real_type y[]
  );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  intersectCircleCircle(
    real_type x1,
    real_type y1,
    real_type theta1,
    real_type kappa1,
    real_type x2,
    real_type y2,
    real_type theta2,
    real_type kappa2,
    real_type s1[],
    real_type s2[]
  );

  /*\
   |   ____        _           ____       ____
   |  / ___|  ___ | |_   _____|___ \__  _|___ \
   |  \___ \ / _ \| \ \ / / _ \ __) \ \/ / __) |
   |   ___) | (_) | |\ V /  __// __/ >  < / __/
   |  |____/ \___/|_| \_/ \___|_____/_/\_\_____|
  \*/
  //!
  //! Class that solve a 2x2 linear system using Pseudo inverse
  //! to manage singular and near singular cases
  //!
  class Solve2x2 {
    integer   i[2], j[2];
    real_type LU[2][2];
    real_type epsi;
    bool      singular;

  public:

    Solve2x2() : epsi(1e-10) {}
    bool factorize( real_type A[2][2] );
    bool solve( real_type const b[2], real_type x[2] ) const;
  };

  //!
  //!  Return the orientation of a triangle
  //!
  //!  \param[in] P1 first point of the triangle
  //!  \param[in] P2 second point of the triangle
  //!  \param[in] P3 third point of the triangle
  //!  \return sign of rotation
  //!
  //!  return +1 = CounterClockwise
  //!  return -1 = Clockwise
  //!  return  0 = flat
  //!
  //!  - CounterClockwise:
  //!    the path P1->P2->P3 turns Counter-Clockwise, i.e.,
  //!    the point P3 is located "on the left" of the line P1-P2.
  //!  - Clockwise:
  //!    the path turns Clockwise, i.e.,
  //!    the point P3 lies "on the right" of the line P1-P2.
  //!  - flat:
  //!    the point P3 is located on the line segment [P1 P2].
  //!
  //!  Algorithm from FileExchage geom2d adapated from Sedgewick's book.
  //!
  integer
  is_counter_clockwise(
    real_type const P1[],
    real_type const P2[],
    real_type const P3[]
  );

  //!
  //! Check if a point is inside a triangle
  //!
  //! \param[in] point point to check if is inside the triangle
  //! \param[in] P1    first point of the triangle
  //! \param[in] P2    second point of the triangle
  //! \param[in] P3    third point of the triangle
  //! \return {0,+1,-1}
  //!         return +1 = Inside
  //!         return -1 = Outsize
  //!         return  0 = on border
  //!

  integer
  is_point_in_triangle(
    real_type const point[],
    real_type const P1[],
    real_type const P2[],
    real_type const P3[]
  );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /*\
   |  __  ____   __  _          _____ _          _
   |  \ \/ /\ \ / / | |_ ___   |_   _| |__   ___| |_ __ _
   |   \  /  \ V /  | __/ _ \    | | | '_ \ / _ \ __/ _` |
   |   /  \   | |   | || (_) |   | | | | | |  __/ || (_| |
   |  /_/\_\  |_|    \__\___/    |_| |_| |_|\___|\__\__,_|
  \*/
  void
  xy_to_guess_angle(
    integer         npts,
    real_type const x[],
    real_type const y[],
    real_type       theta[],
    real_type       theta_min[],
    real_type       theta_max[],
    real_type       omega[],
    real_type       len[]
  );

}

///
/// eof: G2lib.hh
///
