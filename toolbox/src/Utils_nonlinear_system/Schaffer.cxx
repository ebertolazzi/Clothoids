/*\
 |
 |  Author:
 |    Enrico Bertolazzi
 |    University of Trento
 |    Department of Industrial Engineering
 |    Via Sommarive 9, I-38123, Povo, Trento, Italy
 |    email: enrico.bertolazzi@unitn.it
\*/

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

#define SHAFFER_BIBTEX                                                 \
  "@book{brent2013,\n"                                                 \
  "  author    = {Brent, R.P.},\n"                                     \
  "  title     = {Algorithms for Minimization Without Derivatives},\n" \
  "  isbn      = {9780486143682},\n"                                   \
  "  series    = {Dover Books on Mathematics},\n"                      \
  "  year      = {2013},\n"                                            \
  "  publisher = {Dover Publications}\n"                               \
  "}\n"

class SchafferF6 : public NonlinearSystem
{
public:
  SchafferF6() : NonlinearSystem( "Schaffer Function F6", SHAFFER_BIBTEX, 2 ) {}

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type x1 = x( 0 );
    real_type x2 = x( 1 );

    real_type r2 = x1 * x1 + x2 * x2;
    real_type r  = sqrt( r2 );

    if ( r == 0 )
    {
      f( 0 ) = f( 1 ) = 0;
      return;
    }

    // unit vector
    real_type rx1 = x1 / r;
    real_type rx2 = x2 / r;

    real_type D  = 1 + 0.001 * r2;
    real_type D2 = D * D;
    real_type D3 = D2 * D;

    real_type a  = 1 / D2;
    real_type ar = -0.004 * r / D3;

    real_type s  = sin( r );
    real_type b  = s * s - 0.5;
    real_type br = sin( 2 * r );

    real_type S = ar * b + a * br;

    f( 0 ) = S * rx1;
    f( 1 ) = S * rx2;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    real_type x1 = x( 0 );
    real_type x2 = x( 1 );

    real_type r2 = x1 * x1 + x2 * x2;
    real_type r  = sqrt( r2 );

    J.resize( n, n );
    J.setZero();

    if ( r == 0 )
    {
      // corretto: per r=0 il gradiente è 0
      J.makeCompressed();
      return;
    }

    real_type r3 = r2 * r;

    // unit vector
    real_type rx1 = x1 / r;
    real_type rx2 = x2 / r;

    // derivatives of unit vector
    real_type rx1x1 = x2 * x2 / r3;
    real_type rx1x2 = -x1 * x2 / r3;
    real_type rx2x1 = -x1 * x2 / r3;
    real_type rx2x2 = x1 * x1 / r3;

    // auxiliary terms
    real_type D  = 1 + 0.001 * r2;
    real_type D2 = D * D;
    real_type D3 = D2 * D;
    real_type D4 = D2 * D2;

    real_type a  = 1 / D2;
    real_type ar = -0.004 * r / D3;

    // CORRETTO: derivata seconda di a(r)
    real_type arr = -0.004 / D3 + 0.000024 * r2 / D4;

    real_type s   = sin( r );
    real_type b   = s * s - 0.5;
    real_type br  = sin( 2 * r );
    real_type brr = 2 * cos( 2 * r );

    // S = ar*b + a*br
    // dS/dr = arr*b + 2*ar*br + a*brr
    real_type S  = ar * b + a * br;
    real_type Sr = arr * b + 2 * ar * br + a * brr;

    // Jacobiano:
    // J = Sr * (r̂ r̂^T) + S * (d r̂/dx)
    J.insert( 0, 0 ) = Sr * rx1 * rx1 + S * rx1x1;
    J.insert( 0, 1 ) = Sr * rx1 * rx2 + S * rx1x2;
    J.insert( 1, 0 ) = Sr * rx2 * rx1 + S * rx2x1;
    J.insert( 1, 1 ) = Sr * rx2 * rx2 + S * rx2x2;

    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.setZero();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << -5, 10;
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class SchafferF7 : public NonlinearSystem
{
public:
  SchafferF7() : NonlinearSystem( "Schaffer Function F7", SHAFFER_BIBTEX, 2 ) {}

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type x1 = x( 0 );
    real_type x2 = x( 1 );

    real_type r = hypot( x1, x2 );

    if ( r == 0.0 )
    {
      f( 0 ) = f( 1 ) = 0;
      return;
    }

    real_type a  = sqrt( r );
    real_type ar = 0.5 / sqrt( r );

    real_type b  = 1.0 + power2( sin( 50.0 * pow( r, 0.2 ) ) );
    real_type br = 10.0 * sin( 100.0 * pow( r, 0.2 ) ) * pow( r, -0.8 );

    real_type rx1 = x1 / r;
    real_type rx2 = x2 / r;

    f( 0 ) = ( ar * b + a * br ) * rx1;
    f( 1 ) = ( ar * b + a * br ) * rx2;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    real_type x1 = x( 0 );
    real_type x2 = x( 1 );

    real_type r = hypot( x1, x2 );

    if ( r == 0.0 )
    {
      J.setZero();
      return;
    }

    real_type rx1 = x1 / r;
    real_type rx2 = x2 / r;

    real_type r3 = r * r * r;

    real_type rx1x1 = x2 * x2 / r3;
    real_type rx1x2 = -x1 * x2 / r3;
    real_type rx2x1 = -x1 * x2 / r3;
    real_type rx2x2 = x1 * x1 / r3;

    //  F = A * B
    //  dFdX1 = ( Ar * B + A * Br ) * Rx1
    //  d2FdX1dX1 = ( Arr * B + Ar * Br ) * Rx1^2 + ( Ar * B + A * Br ) * Rx1x1
    //  etc
    real_type a   = sqrt( r );
    real_type ar  = 0.5 / sqrt( r );
    real_type arr = -0.25 / sqrt( r3 );

    real_type b   = 1.0 + power2( sin( 50.0 * pow( r, 0.2 ) ) );
    real_type br  = 10.0 * sin( 100.0 * pow( r, 0.2 ) ) * pow( r, -0.8 );
    real_type brr = 200.0 * cos( 100.0 * pow( r, 0.2 ) ) * pow( r, -1.6 ) -
                    10.0 * sin( 100.0 * pow( r, 0.2 ) ) * 0.8 * pow( r, -1.8 );

    J.resize( n, n );
    J.setZero();

    J.insert( 0, 0 ) = ( arr * b + 2 * ar * br + a * brr ) * rx1 * rx1 + ( ar * b + a * br ) * rx1x1;
    J.insert( 0, 1 ) = ( arr * b + 2 * ar * br + a * brr ) * rx1 * rx2 + ( ar * b + a * br ) * rx1x2;
    J.insert( 1, 0 ) = ( arr * b + 2 * ar * br + a * brr ) * rx2 * rx1 + ( ar * b + a * br ) * rx2x1;
    J.insert( 1, 1 ) = ( arr * b + 2 * ar * br + a * brr ) * rx2 * rx2 + ( ar * b + a * br ) * rx2x2;

    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.setZero();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << -5, 10;
  }
};
