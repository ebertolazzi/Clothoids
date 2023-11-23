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
/// file: Clothoid_compatibility.hxx
///

//!
//! Clothoid curve total variation of the angle.
//!
real_type
thetaTotalVariation() const
{ return theta_total_variation(); }

//!
//! Max and min angle of the curve.
//!
real_type
thetaMinMax( real_type & thMin, real_type & thMax ) const
{ return theta_min_max( thMin, thMax ); }

//!
//! Clothoid angle range.
//!
real_type
deltaTheta() const
{ return deltaTheta(); }

//!
//! Max and min of the curvatire of the clothoid curve.
//!
real_type
curvatureMinMax( real_type & kMin, real_type & kMax ) const
{ return curvature_min_max( kMin, kMax ); }

//!
//! Clothoid total curvature variation.
//!
real_type
curvatureTotalVariation() const
{ return curvature_total_variation(); }

//!
//! Given the clothoid curve \f$ P(s) \f$ compute.
//!
//! \f[
//!    \int_0^L |P''(s)|^2 \mathrm{d}s
//! \f]
//!
real_type
integralCurvature2() const
{ return integral_curvature2(); }

//!
//! Given the clothoid curve \f$ P(s) \f$ compute.
//!
//! \f[
//!    \int_0^L |P'''(s)|^2 \mathrm{d}s
//! \f]
//!
real_type
integralJerk2() const
{ return integral_jerk2(); }

//!
//! Given the clothoid curve \f$ P(s) \f$ compute.
//!
//! \f[
//!    \int_0^L |P''''(s)|^2 \mathrm{d}s
//! \f]
//!
real_type
integralSnap2() const
{ return integral_snap2(); }

real_type
distanceBySample(
  real_type   ds,
  real_type   qx,
  real_type   qy,
  real_type & S
) const {
  real_type X, Y;
  return closest_point_by_sample( ds, qx, qy, X, Y, S );
}

real_type
distanceBySample(
  real_type ds,
  real_type qx,
  real_type qy
) const {
  real_type X, Y, S;
  return closest_point_by_sample( ds, qx, qy, X, Y, S );
}

void
evaluate(
  real_type   s,
  real_type & th,
  real_type & k,
  real_type & x,
  real_type & y
) const override
{ m_CD.evaluate( s, th, k, x, y ); }

void
changeOrigin( real_type newx0, real_type newy0 )
{ change_origin( newx0, newy0 ); }

void
changeCurvilinearOrigin( real_type s0, real_type newL )
{ change_curvilinear_origin( s0, newL ); }

real_type thetaBegin()                 const { return theta_begin(); }
real_type thetaEnd()                   const { return theta_end(); }
real_type kappaBegin()                 const { return kappa_begin(); }
real_type kappaEnd()                   const { return kappa_end(); }
real_type xBegin()                     const { return x_begin(); }
real_type yBegin()                     const { return y_begin(); }
real_type xEnd()                       const { return x_end(); }
real_type yEnd()                       const { return y_end(); }
real_type xBegin_ISO( real_type offs ) const { return x_begin_ISO( offs ); }
real_type yBegin_ISO( real_type offs ) const { return y_Begin_ISO( offs ); }
real_type xEnd_ISO( real_type offs )   const { return x_end_ISO( offs ); }
real_type yEnd_ISO( real_type offs )   const { return y_end_ISO( offs ); }

real_type
closestPointBySample(
  real_type   ds,
  real_type   qx,
  real_type   qy,
  real_type & X,
  real_type & Y,
  real_type & S
) const {
  return closest_point_by_sample( ds, qx, qy, X, Y, S );
}

integer
closestPoint_ISO(
  real_type   qx,
  real_type   qy,
  real_type & x,
  real_type & y,
  real_type & s,
  real_type & t,
  real_type & dst
) const {
  return closest_point_ISO( qx, qy, x, y, s, t, dst );
}

integer
closestPoint_ISO(
  real_type   qx,
  real_type   qy,
  real_type   offs,
  real_type & x,
  real_type & y,
  real_type & s,
  real_type & t,
  real_type & dst
) const {
  return closest_point_ISO( qx, qy, offs, x, y, s, t, dst );
}

real_type
distanceBySample(
  real_type   ds,
  real_type   qx,
  real_type   qy,
  real_type & S
) const {
  real_type X, Y;
  return closest_point_by_sample( ds, qx, qy, X, Y, S );
}

real_type
distanceBySample(
  real_type ds,
  real_type qx,
  real_type qy
) const {
  real_type X, Y, S;
  return closest_point_by_sample( ds, qx, qy, X, Y, S );
}

///
/// eof: Clothoid_compatibulity.hxx
///
