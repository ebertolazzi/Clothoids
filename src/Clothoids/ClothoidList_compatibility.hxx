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

real_type tx_Begin() const { return tx_begin(); }
real_type ty_Begin() const { return ty_begin(); }
real_type tx_End()   const { return tx_end(); }
real_type ty_End()   const { return ty_end(); }

real_type nx_Begin_ISO() const { return nx_begin_ISO(); }
real_type ny_Begin_ISO() const { return ny_begin_ISO(); }
real_type nx_End_ISO()   const { return nx_end_ISO(); }
real_type ny_End_ISO()   const { return ny_end_ISO(); }

//!
//! Return the clothoid list as a list of nodes and curvatures
//!
//! \param[out] s     nodes
//! \param[out] kappa curvature
//!
void
getSK( real_type * s, real_type * kappa ) const
{ get_SK( s, kappa ); }

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
{ get_SK( s, kappa ); }

//!
//! Return the clothoid list as a list of nodes angles and curvatures
//!
//! \param[out] s     nodes
//! \param[out] theta angles
//! \param[out] kappa curvature
//!
void
getSTK(
  real_type * s,
  real_type * theta,
  real_type * kappa
) const
{ get_STK( s, theta, kappa ); }

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
{ get_STK( s, theta, kappa ); }

//!
//! Return the points of the clothoid list at breakpoints
//!
//! \param[out] x x-coordinates
//! \param[out] y y-coordinates
//!
void
getXY( real_type * x, real_type * y ) const
{ get_XY( x, y ); }

void
getDeltaTheta( real_type * delta_theta ) const
{ get_delta_theta( delta_theta ); }

void
getDeltaKappa( real_type * deltaKappa ) const
{ get_delta_kappa( deltaKappa ); }

void
changeOrigin( real_type newx0, real_type newy0 )
{ change_origin( newx0, newy0 ); }

real_type thetaBegin()                 const { return theta_begin(); }
real_type thetaEnd()                   const { return theta_end(); }
real_type xBegin()                     const { return x_begin(); }
real_type yBegin()                     const { return y_begin(); }
real_type xEnd()                       const { return x_end(); }
real_type yEnd()                       const { return y_end(); }
real_type xBegin_ISO( real_type offs ) const { return x_begin_ISO( offs ); }
real_type yBegin_ISO( real_type offs ) const { return y_begin_ISO( offs ); }
real_type xEnd_ISO( real_type offs )   const { return x_end_ISO( offs ); }
real_type yEnd_ISO( real_type offs )   const { return y_end_ISO( offs ); }

integer numSegments() const { return num_segments(); }

integer
closestSegment( real_type qx, real_type qy ) const {
  return closest_segment( qx, qy );
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
  return closest_point_in_range_ISO(
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
  return closest_point_in_range_SAE(
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
  return closest_point_in_s_range_ISO(
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
  return closest_point_in_s_range_SAE(
    qx, qy, s_begin, s_end, x, y, s, t, dst, icurve
  );
}

///
/// eof: ClothoidList_compatibility.hxx
///
