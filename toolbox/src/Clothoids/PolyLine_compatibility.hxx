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
/// file: PolyLine_compatibility.hxx
///

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
closestPoint_ISO(
  real_type   x,
  real_type   y,
  real_type & X,
  real_type & Y,
  real_type & S,
  real_type & T,
  real_type & DST
) const {
  return this->closest_point_ISO( x, y, X, Y, S, T, DST );
}

///
/// eof: PolyLine_compatibility.hxx
///
