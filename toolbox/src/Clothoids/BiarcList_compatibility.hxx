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
/// file: BiarcList_compatibility.hxx
///

void
getSTK(
  real_type s[],
  real_type theta[],
  real_type kappa[]
) const
{ this->get_STK( s, theta, kappa ); }

//!
//! Return the biarc XY nodes
//!
//! \param[out] x \f$x\f$-nodes
//! \param[out] y \f$y\f$-nodes
//
void
getXY( real_type x[], real_type y[] ) const
{ this->get_XY( x, y ); }

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

integer numSegments() const { return this->num_segments(); }

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

///
/// eof: BiarcList_compatibility.hxx
///
