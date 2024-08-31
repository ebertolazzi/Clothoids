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
/// file: ClothoidList_compatibility.hxx
///

//real_type curvatureTotalVariation() const { return this->curvature_total_variation(); }
//real_type integralCurvature2() const { return this->integral_curvature2(); }
//real_type integralJerk2() const { return this->integral_jerk2(); }
//real_type integralSnap2() const { return this->integral_snap2(); }

real_type tx_Begin() const { return this->tx_begin(); }
real_type ty_Begin() const { return this->ty_begin(); }
real_type tx_End()   const { return this->tx_end(); }
real_type ty_End()   const { return this->ty_end(); }

real_type nx_Begin_ISO() const { return this->nx_begin_ISO(); }
real_type ny_Begin_ISO() const { return this->ny_begin_ISO(); }
real_type nx_End_ISO()   const { return this->nx_end_ISO(); }
real_type ny_End_ISO()   const { return this->ny_end_ISO(); }

//!
//! Return the clothoid list as a list of nodes and curvatures
//!
//! \param[out] s     nodes
//! \param[out] kappa curvature
//!
void
getSK( real_type s[], real_type kappa[] ) const
{ this->get_SK( s, kappa ); }

//!
//! Return the clothoid list as a list of nodes and curvatures
//!
//! \param[out] s     nodes
//! \param[out] kappa curvature
//!
void
getSK(
  vector<real_type> & s,
  vector<real_type> & kappa
) const
{ this->get_SK( s, kappa ); }

//!
//! Return the clothoid list as a list of nodes angles and curvatures
//!
//! \param[out] s     nodes
//! \param[out] theta angles
//! \param[out] kappa curvature
//!
void
getSTK(
  real_type s[],
  real_type theta[],
  real_type kappa[]
) const
{ this->get_STK( s, theta, kappa ); }

//!
//! Return the clothoid list as a list of nodes angles and curvatures
//!
//! \param[out] s     nodes
//! \param[out] theta angles
//! \param[out] kappa curvature
//!
void
getSTK(
  vector<real_type> & s,
  vector<real_type> & theta,
  vector<real_type> & kappa
) const
{ this->get_STK( s, theta, kappa ); }

//!
//! Return the points of the clothoid list at breakpoints
//!
//! \param[out] x \f$x\f$-coordinates
//! \param[out] y \f$y\f$-coordinates
//!
void
getXY( real_type x[], real_type y[] ) const
{ this->get_XY( x, y ); }

void
getDeltaTheta( real_type delta_theta[] ) const
{ this->get_delta_theta( delta_theta ); }

void
getDeltaKappa( real_type deltaKappa[] ) const
{ this->get_delta_kappa( deltaKappa ); }

void
changeOrigin( real_type newx0, real_type newy0 )
{ this->change_origin( newx0, newy0 ); }

real_type thetaBegin()                 const { return this->theta_begin(); }
real_type thetaEnd()                   const { return this->theta_end(); }
real_type xBegin()                     const { return this->x_begin(); }
real_type yBegin()                     const { return this->y_begin(); }
real_type xEnd()                       const { return this->x_end(); }
real_type yEnd()                       const { return this->y_end(); }
real_type xBegin_ISO( real_type offs ) const { return this->x_begin_ISO( offs ); }
real_type yBegin_ISO( real_type offs ) const { return this->y_begin_ISO( offs ); }
real_type xEnd_ISO( real_type offs )   const { return this->x_end_ISO( offs ); }
real_type yEnd_ISO( real_type offs )   const { return this->y_end_ISO( offs ); }

integer numSegments() const { return num_segments(); }

integer
closestSegment( real_type qx, real_type qy ) const {
  return this->closest_segment( qx, qy );
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
  return this->closest_point_ISO( qx, qy, x, y, s, t, dst );
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
  return this->closest_point_ISO( qx, qy, offs, x, y, s, t, dst );
}

integer
closestPointInRange_ISO(
  real_type   qx,
  real_type   qy,
  integer     icurve_begin,
  integer     icurve_end,
  real_type & x,
  real_type & y,
  real_type & s,
  real_type & t,
  real_type & dst,
  integer   & icurve
) const {
  return this->closest_point_in_range_ISO(
    qx, qy, icurve_begin, icurve_end, x, y, s, t, dst, icurve
  );
}
integer
closestPointInRange_SAE(
  real_type   qx,
  real_type   qy,
  integer     icurve_begin,
  integer     icurve_end,
  real_type & x,
  real_type & y,
  real_type & s,
  real_type & t,
  real_type & dst,
  integer   & icurve
) const {
  return this->closest_point_in_range_SAE(
    qx, qy, icurve_begin, icurve_end, x, y, s, t, dst, icurve
  );
}

integer
closestPointInSRange_ISO(
  real_type   qx,
  real_type   qy,
  real_type   s_begin,
  real_type   s_end,
  real_type & x,
  real_type & y,
  real_type & s,
  real_type & t,
  real_type & dst,
  integer   & icurve
) const {
  return this->closest_point_in_s_range_ISO(
    qx, qy, s_begin, s_end, x, y, s, t, dst, icurve
  );
}

integer
closestPointInSRange_SAE(
  real_type   qx,
  real_type   qy,
  real_type   s_begin,
  real_type   s_end,
  real_type & x,
  real_type & y,
  real_type & s,
  real_type & t,
  real_type & dst,
  integer   & icurve
) const {
  return this->closest_point_in_s_range_SAE(
    qx, qy, s_begin, s_end, x, y, s, t, dst, icurve
  );
}

///
/// eof: ClothoidList_compatibility.hxx
///
