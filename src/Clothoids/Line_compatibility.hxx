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
/// file: Line_compatibility.hxx
///

void
changeOrigin( real_type newx0, real_type newy0 )
{ change_origin( newx0, newy0 ); }

real_type xBegin()                     const { return x_begin(); }
real_type yBegin()                     const { return y_begin(); }
real_type xEnd()                       const { return x_end(); }
real_type yEnd()                       const { return y_end(); }
real_type xBegin_ISO( real_type offs ) const { return x_begin_ISO( offs ); }
real_type yBegin_ISO( real_type offs ) const { return y_begin_ISO( offs ); }
real_type xEnd_ISO( real_type offs )   const { return x_end_ISO( offs ); }
real_type yEnd_ISO( real_type offs )   const { return y_end_ISO( offs ); }

real_type tx_Begin()     const { return m_c0; }
real_type ty_Begin()     const { return m_s0; }
real_type tx_End()       const { return m_c0; }
real_type ty_End()       const { return m_s0; }
real_type nx_Begin_ISO() const { return -m_s0; }
real_type ny_Begin_ISO() const { return m_c0; }
real_type nx_End_ISO()   const { return -m_s0; }
real_type ny_End_ISO()   const { return m_c0; }

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

///
/// eof: Line_compatibility.hxx
///
