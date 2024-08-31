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
/// file: BaseCurve_compatibility.hh
///

//!
//! Length of the curve with offset (ISO/SAE)
//!
real_type
length( real_type offs ) const
{ return G2lib::use_ISO ? this->length_ISO(offs) : this->length_SAE(offs); }

//!
//! Compute the bounding box of the curve (ISO/SAE).
//!
//! \param[in]  offs curve offset
//! \param[out] xmin left bottom
//! \param[out] ymin left bottom
//! \param[out] xmax right top
//! \param[out] ymax right top
//!
void
bbox(
  real_type   offs,
  real_type & xmin,
  real_type & ymin,
  real_type & xmax,
  real_type & ymax
) const {
  if ( G2lib::use_ISO ) this->bbox_ISO( offs, xmin, ymin, xmax, ymax );
  else                  this->bbox_SAE( offs, xmin, ymin, xmax, ymax );
}

//!
//! Initial \f$x\f$-coordinate with offset (ISO or SAE).
//!
real_type xBegin( real_type offs ) const
{ return G2lib::use_ISO ? this->x_begin_ISO(offs) : this->x_begin_SAE(offs); }

//!
//! Initial \f$y\f$-coordinate with offset (ISO or SAE).
//!
real_type yBegin( real_type offs ) const
{ return G2lib::use_ISO ? this->y_begin_ISO(offs) : this->y_begin_SAE(offs); }

//!
//! Final \f$x\f$-coordinate with offset (ISO or SAE).
//!
real_type xEnd( real_type offs ) const
{ return G2lib::use_ISO ? this->x_end_ISO(offs) : this->x_end_SAE(offs); }

//!
//! Final \f$y\f$-coordinate with offset (ISO or SAE).
//!
real_type yEnd( real_type offs ) const
{ return G2lib::use_ISO ? this->y_end_ISO(offs) : this->y_end_SAE(offs); }

//!
//! Intial normal \f$x\f$-coordinate.
//!
real_type nx_Begin() const
{ return G2lib::use_ISO ? this->nx_begin_ISO() : this->nx_begin_SAE(); }

//!
//! Intial normal \f$y\f$-coordinate.
//!
real_type ny_Begin() const
{ return G2lib::use_ISO ? this->ny_begin_ISO() : this->ny_begin_SAE(); }

//!
//! Final normal \f$x\f$-coordinate.
//!
real_type nx_End() const
{ return G2lib::use_ISO ? this->nx_end_ISO() : this->nx_end_SAE(); }

//!
//! Final normal \f$y\f$-coordinate.
//!
real_type ny_End() const
{ return G2lib::use_ISO ? this->ny_end_ISO() : this->ny_end_SAE(); }

//!
//! Normal \f$x\f$-coordinate at curvilinear coordinate \f$s\f$ (ISO/SAE).
//!
real_type nx( real_type s ) const
{ return G2lib::use_ISO ? this->nx_ISO(s) : this->nx_SAE(s); }

//!
//! Normal derivative \f$x\f$-coordinate at curvilinear coordinate \f$s\f$ (ISO/SAE).
//!
real_type nx_D( real_type s ) const
{ return G2lib::use_ISO ? this->nx_ISO_D(s) : this->nx_SAE_D(s); }

//!
//! Normal second derivative \f$x\f$-coordinate at curvilinear coordinate \f$s\f$ (ISO/SAE).
//!
real_type nx_DD( real_type s ) const
{ return G2lib::use_ISO ? this->nx_ISO_DD(s) : this->nx_SAE_DD(s); }

//!
//! Normal third derivative \f$x\f$-coordinate at curvilinear coordinate \f$s\f$ (ISO/SAE).
//!
real_type nx_DDD( real_type s ) const
{ return G2lib::use_ISO ? this->nx_ISO_DDD(s) : this->nx_SAE_DDD(s); }

//!
//! Normal \f$y\f$-coordinate at curvilinear coordinate \f$s\f$ (ISO/SAE).
//!
real_type ny( real_type s ) const
{ return G2lib::use_ISO ? this->ny_ISO(s) : this->ny_SAE(s); }

//!
//! Normal derivative \f$y\f$-coordinate at curvilinear coordinate \f$s\f$ (ISO/SAE).
//!
real_type ny_D( real_type s ) const
{ return G2lib::use_ISO ? this->ny_ISO_D(s) : this->ny_SAE_D(s); }

//!
//! Normal second derivative \f$y\f$-coordinate at curvilinear coordinate \f$s\f$ (ISO/SAE).
//!
real_type ny_DD( real_type s ) const
{ return G2lib::use_ISO ? this->ny_ISO_DD(s) : this->ny_SAE_DD(s); }

//!
//! Normal third derivative \f$y\f$-coordinate at curvilinear coordinate \f$s\f$ (ISO/SAE).
//!
real_type ny_DDD( real_type s ) const
{ return G2lib::use_ISO ? this->ny_ISO_DDD(s) : this->ny_SAE_DDD(s); }

//!
//! Initial tangent \f$x\f$-coordinate.
//!
real_type tx_Begin() const { return this->tx_begin(); }

//!
//! Initial tangent \f$y\f$-coordinate.
//!
real_type ty_Begin() const { return this->ty_begin(); }

//!
//! Final tangent \f$x\f$-coordinate.
//!
real_type tx_End() const { return this->tx_end(); }

//!
//! Final tangent \f$y\f$-coordinate.
//!
real_type ty_End() const { return this->ty_end(); }

//!
//! Intial normal \f$x\f$-coordinate (ISO).
//!
real_type nx_Begin_ISO() const { return nx_begin_ISO(); }

//!
//! Intial normal \f$y\f$-coordinate (ISO).
//!
real_type ny_Begin_ISO() const { return this->ny_begin_ISO(); }

//!
//! Final normal \f$x\f$-coordinate (ISO).
//!
real_type nx_End_ISO() const { return nx_end_ISO(); }

//!
//! Final normal \f$y\f$-coordinate (ISO).
//!
real_type ny_End_ISO() const { return this->ny_end_ISO(); }

//!
//! Intial normal \f$x\f$-coordinate (SAE).
//!
real_type nx_Begin_SAE() const { return nx_begin_SAE(); }

//!
//! Intial normal \f$y\f$-coordinate (SAE).
//!
real_type ny_Begin_SAE() const { return ny_begin_SAE(); }

//!
//! Final normal \f$x\f$-coordinate (SAE).
//!
real_type nx_End_SAE() const { return nx_end_SAE(); }

//!
//! Intial normal \f$y\f$-coordinate (SAE).
//!
real_type ny_End_SAE() const { return ny_end_SAE(); }

//!
//! Normal at curvilinear coordinate \f$s\f$ (ISO/SAE).
//!
void
nor( real_type s, real_type & nx, real_type & ny ) const {
  if ( G2lib::use_ISO ) this->nor_ISO(s,nx,ny);
  else                  this->nor_SAE(s,nx,ny);
}

//!
//! Normal derivative at curvilinear coordinate \f$s\f$ (ISO/SAE).
//!
void
nor_D( real_type s, real_type & nx_D, real_type & ny_D ) const {
  if ( G2lib::use_ISO ) this->nor_ISO_D(s,nx_D,ny_D);
  else                  this->nor_SAE_D(s,nx_D,ny_D);
}

//!
//! Normal second derivative at curvilinear coordinate \f$s\f$ (ISO/SAE).
//!
void
nor_DD( real_type s, real_type & nx_DD, real_type & ny_DD ) const {
  if ( G2lib::use_ISO ) this->nor_ISO_DD(s,nx_DD,ny_DD);
  else                  this->nor_SAE_DD(s,nx_DD,ny_DD);
}

//!
//! Normal third at curvilinear coordinate \f$s\f$ (ISO/SAE).
//!
void
nor_DDD( real_type s, real_type & nx_DDD, real_type & ny_DDD ) const {
  if ( G2lib::use_ISO ) this->nor_ISO_DDD(s,nx_DDD,ny_DDD);
  else                  this->nor_SAE_DDD(s,nx_DDD,ny_DDD);
}

//!
//! Evaluate curve with offset at curvilinear coordinate \f$s\f$ (ISO/SAE).
//!
//! \param[in]  s    curvilinear coordinate
//! \param[in]  offs offset
//! \param[out] th   angle
//! \param[out] k    curvature
//! \param[out] x    \f$x\f$-coordinate
//! \param[out] y    \f$y\f$-coordinate
//!
void
evaluate(
  real_type   s,
  real_type   offs,
  real_type & th,
  real_type & k,
  real_type & x,
  real_type & y
) const {
  if ( G2lib::use_ISO ) this->evaluate_ISO( s, offs, th, k, x, y );
  else                  this->evaluate_SAE( s, offs, th, k, x, y );
}

//!
//! \f$x\f$-coordinate at curvilinear coordinate \f$s\f$ with offset `offs` (ISO/SAE).
//!
real_type
X( real_type s, real_type offs ) const
{ return G2lib::use_ISO ? this->X_ISO( s, offs ) : this->X_SAE( s, offs ); }

//!
//! \f$y\f$-coordinate at curvilinear coordinate \f$s\f$ with offset `offs` (ISO/SAE).
//!
real_type
Y( real_type s, real_type offs ) const
{ return G2lib::use_ISO ? this->Y_ISO( s, offs ) : this->Y_SAE( s, offs ); }

//!
//! \f$x\f$-coordinate derivative at curvilinear coordinate \f$s\f$ with offset `offs` (ISO/SAE).
//!
real_type
X_D( real_type s, real_type offs ) const
{ return G2lib::use_ISO ? this->X_ISO_D( s, offs ) : this->X_SAE_D( s, offs ); }

//!
//! \f$y\f$-coordinate derivative at curvilinear coordinate \f$s\f$ with offset `offs` (ISO/SAE).
//!
real_type
Y_D( real_type s, real_type offs ) const
{ return G2lib::use_ISO ? this->Y_ISO_D( s, offs ) : this->Y_SAE_D( s, offs ); }

//!
//! \f$x\f$-coordinate second derivative at curvilinear coordinate \f$s\f$ with offset `offs` (ISO/SAE).
//!
real_type
X_DD( real_type s, real_type offs ) const
{ return G2lib::use_ISO ? this->X_ISO_DD( s, offs ) : this->X_SAE_DD( s, offs ); }

//!
//! \f$y\f$-coordinate second derivative at curvilinear coordinate \f$s\f$ with offset `offs` (ISO/SAE).
//!
real_type
Y_DD( real_type s, real_type offs ) const
{ return G2lib::use_ISO ? this->Y_ISO_DD( s, offs ) : this->Y_SAE_DD( s, offs ); }

//!
//! \f$x\f$-coordinate third derivative at curvilinear coordinate \f$s\f$ with offset `offs` (ISO/SAE).
//!
real_type
X_DDD( real_type s, real_type offs ) const
{ return G2lib::use_ISO ? this->X_ISO_DDD( s, offs ) : this->X_SAE_DDD( s, offs ); }

//!
//! \f$y\f$-coordinate third derivative at curvilinear coordinate \f$s\f$ with offset `offs` (ISO/SAE).
//!
real_type
Y_DDD( real_type s, real_type offs ) const
{ return G2lib::use_ISO ? this->Y_ISO_DDD( s, offs ) : this->Y_SAE_DDD( s, offs ); }

//!
//! Compute curve at position `s` with offset `offs` (ISO/SAE).
//!
//! \param[in]  s     parameter on the curve
//! \param[in]  offs  offset of the curve
//! \param[out] x     \f$ x\f$-coordinate
//! \param[out] y     \f$ y\f$-coordinate
//!
void
eval(
  real_type   s,
  real_type   offs,
  real_type & x,
  real_type & y
) const {
  if ( G2lib::use_ISO ) this->eval_ISO( s, offs, x, y );
  else                  this->eval_SAE( s, offs, x, y );
}

//!
//! Compute derivative curve at position `s` with offset `offs`  (ISO/SAE).
//!
//! \param[in]  s     parameter on the curve
//! \param[in]  offs  offset of the curve
//! \param[out] x_D   \f$x\f$-coordinate first derivative
//! \param[out] y_D   \f$y\f$-coordinate first derivative
//!
void
eval_D(
  real_type   s,
  real_type   offs,
  real_type & x_D,
  real_type & y_D
) const {
  if ( G2lib::use_ISO ) this->eval_ISO_D( s, offs, x_D, y_D );
  else                  this->eval_SAE_D( s, offs, x_D, y_D );
}

//!
//! Compute second derivative curve at position `s` with offset `offs`  (ISO/SAE).
//!
//! \param[in]  s     parameter on the curve
//! \param[in]  offs  offset of the curve
//! \param[out] x_DD  \f$x\f$-coordinate second derivative
//! \param[out] y_DD  \f$y\f$-coordinate second derivative
//!
void
eval_DD(
  real_type   s,
  real_type   offs,
  real_type & x_DD,
  real_type & y_DD
) const {
  if ( G2lib::use_ISO ) this->eval_ISO_DD( s, offs, x_DD, y_DD );
  else                  this->eval_SAE_DD( s, offs, x_DD, y_DD );
}

//!
//! Compute third derivative curve at position `s` with offset `offs` (ISO/SAE).
//!
//! \param[in]  s     parameter on the curve
//! \param[in]  offs  offset of the curve
//! \param[out] x_DDD \f$x\f$-coordinate third derivative
//! \param[out] y_DDD \f$y\f$-coordinate third derivative
//!
void
eval_DDD(
  real_type   s,
  real_type   offs,
  real_type & x_DDD,
  real_type & y_DDD
) const {
  if ( G2lib::use_ISO ) this->eval_ISO_DDD( s, offs, x_DDD, y_DDD );
  else                  this->eval_SAE_DDD( s, offs, x_DDD, y_DDD );
}

void
changeOrigin( real_type newx0, real_type newy0 )
{ change_origin( newx0, newy0 ); }

real_type thetaBegin()                 const { return this->theta_begin(); }
real_type thetaEnd()                   const { return this->theta_end(); }
real_type kappaBegin()                 const { return this->kappa_begin(); }
real_type kappaEnd()                   const { return this->kappa_end(); }
real_type xBegin()                     const { return this->x_begin(); }
real_type yBegin()                     const { return this->y_begin(); }
real_type xEnd()                       const { return this->x_end(); }
real_type yEnd()                       const { return this->y_end(); }
real_type xBegin_ISO( real_type offs ) const { return this->x_begin_ISO(offs); }
real_type yBegin_ISO( real_type offs ) const { return this->y_begin_ISO(offs); }
real_type xEnd_ISO( real_type offs )   const { return this->x_end_ISO(offs); }
real_type yEnd_ISO( real_type offs )   const { return this->y_end_ISO(offs); }
real_type xBegin_SAE( real_type offs ) const { return this->x_begin_SAE(offs); }
real_type yBegin_SAE( real_type offs ) const { return this->y_begin_SAE(offs); }
real_type xEnd_SAE( real_type offs )   const { return this->x_end_SAE(offs); }
real_type yEnd_SAE( real_type offs )   const { return this->y_end_SAE(offs); }

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
closestPoint_SAE(
  real_type   qx,
  real_type   qy,
  real_type & x,
  real_type & y,
  real_type & s,
  real_type & t,
  real_type & dst
) const {
  return this->closest_point_SAE( qx, qy, x, y, s, t, dst );
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
closestPoint_SAE(
  real_type   qx,
  real_type   qy,
  real_type   offs,
  real_type & x,
  real_type & y,
  real_type & s,
  real_type & t,
  real_type & dst
) const {
  return this->closest_point_SAE( qx, qy, offs, x, y, s, t, dst );
}

///
/// eof: BaseCurve_compatibility.hxx
///
