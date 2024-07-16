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
/// file: Circle_compatibility.hxx
///

real_type
thetaTotalVariation() const
{ return this->theta_total_variation(); }

real_type
thetaMinMax( real_type & thMin, real_type & thMax ) const
{ return this->theta_min_max(thMin,thMax); }

void
changeOrigin( real_type newx0, real_type newy0 )
{ this->change_origin( newx0, newy0 ); }

void
changeCurvilinearOrigin( real_type s0, real_type newL )
{ this->change_curvilinear_origin( s0, newL ); }

real_type thetaBegin()                 const { return this->theta_begin(); }
real_type thetaEnd()                   const { return this->theta_end(); }
real_type kappaBegin()                 const { return this->kappa_begin(); }
real_type kappaEnd()                   const { return this->kappa_end(); }
real_type xBegin()                     const { return this->x_begin(); }
real_type yBegin()                     const { return this->y_begin(); }
real_type xEnd()                       const { return this->x_end(); }
real_type yEnd()                       const { return this->y_end(); }
real_type xBegin_ISO( real_type offs ) const { return this->x_begin_ISO( offs ); }
real_type yBegin_ISO( real_type offs ) const { return this->y_begin_ISO( offs ); }
real_type xEnd_ISO( real_type offs )   const { return this->x_end_ISO( offs ); }
real_type yEnd_ISO( real_type offs )   const { return this->y_end_ISO( offs ); }

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
//!
//! Return \f$ \sin \theta_0 \f$ where
//! \f$ \theta_0 \f$ is the initial tangent angle.
//!
real_type sinTheta0() const { return sin(m_theta0); }

//!
//! Return \f$ \cos \theta_0 \f$ where
//! \f$ \theta_0 \f$ is the initial tangent angle.
//!
real_type cosTheta0() const { return cos(m_theta0); }

//!
//! Return the length of the arc that
//! can approximated by a line segment.
//!
real_type lenTolerance( real_type tol ) const
{ return this->len_tolerance( tol ); }

///
/// eof: Circle_compatibility.hxx
///
