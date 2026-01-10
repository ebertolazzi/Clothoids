/*\
 |
 |  Author:
 |    Enrico Bertolazzi
 |    University of Trento
 |    Department of Industrial Engineering
 |    Via Sommarive 9, I-38123, Povo, Trento, Italy
 |    email: enrico.bertolazzi@unitn.it
\*/

#define SINGULAR_SYSTEM_BIBTEX                                                 \
  "@article{Waziri:2011,\n"                                                    \
  "  author  = {Waziri Yusuf, Mohammed and June, Leong Wah and Hassan, Malik " \
  "Abu},\n"                                                                    \
  "  title   = {Jacobian-free diagonal {N}ewton's method for solving "         \
  "nonlinear\n"                                                                \
  "             systems with singular {J}acobian},\n"                          \
  "  joirnal = {Malaysian Journal of Mathematical Sciences},\n"                \
  "  volume  = {5},\n"                                                         \
  "  year    = {2011},\n"                                                      \
  "  number  = {2},\n"                                                         \
  "  pages   = {241--255}\n"                                                   \
  "}\n"

/*
  ￼Modified Newton’s method for systems of nonlinear equations with singular
  Jacobian
*/

// Problem 1. (Jose et al. (2009))

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class SingularSystemA : public NonlinearSystem
{
public:
  SingularSystemA() : NonlinearSystem( "Singular System A", SINGULAR_SYSTEM_BIBTEX, 2 ) {}

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = power2( x( 0 ) - 1 ) * ( x( 0 ) - x( 1 ) );
    f( 1 ) = power5( x( 1 ) - 2 ) * cos( 2 * x( 0 ) / x( 1 ) );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type t      = 2 * x( 0 ) / x( 1 );
    real_type S      = sin( t );
    real_type C      = cos( t );
    J.insert( 0, 0 ) = ( x( 0 ) - 1 ) * ( 3 * x( 0 ) - 2 * x( 1 ) - 1 );
    J.insert( 0, 1 ) = -power2( x( 0 ) - 1 );
    J.insert( 1, 0 ) = -2 * power5( x( 1 ) - 2 ) * S / x( 1 );
    J.insert( 1, 1 ) = power4( x( 1 ) - 2 ) * ( 5 * C + 2 * ( x( 1 ) - 2 ) * x( 0 ) * S / ( x( 1 ) * x( 1 ) ) );
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 1, 2;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 4 );
    auto & x0{ x_vec[0] };
    auto & x1{ x_vec[1] };
    auto & x2{ x_vec[2] };
    auto & x3{ x_vec[3] };
    x0.resize( n );
    x1.resize( n );
    x2.resize( n );
    x3.resize( n );

    x0 << 1.5, 2.5;
    x1 << 2, 5;
    x2 << 0, 3;
    x3 << 0.5, 2;
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class SingularSystemB : public NonlinearSystem
{
public:
  SingularSystemB() : NonlinearSystem( "Singular System B", SINGULAR_SYSTEM_BIBTEX, 3 ) {}

  virtual void evaluate( Vector const & X, Vector & f ) const override
  {
    real_type x = X( 0 );
    real_type y = X( 1 );
    real_type z = X( 2 );
    f( 0 )      = power4( x - 1 ) * exp( y );
    f( 1 )      = power5( y - 2 ) * ( x * y - 1 );
    f( 2 )      = power6( z + 4 );
  }

  virtual void jacobian( Vector const & X, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type x      = X( 0 );
    real_type y      = X( 1 );
    real_type z      = X( 2 );
    J.insert( 0, 0 ) = 4 * power3( x - 1 ) * exp( y );
    J.insert( 0, 1 ) = power4( x - 1 ) * exp( y );
    J.insert( 1, 0 ) = power5( y - 2 ) * y;
    J.insert( 1, 1 ) = power4( y - 2 ) * ( 6 * x * y - 2 * x - 5 );
    J.insert( 2, 2 ) = 6 * power5( z + 4 );
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 1, 2, -4;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 2 );
    auto & x0{ x_vec[0] };
    auto & x1{ x_vec[1] };
    x0.resize( n );
    x1.resize( n );

    x0 << 2, 1, -2;
    x1 << 4, 5, 6;
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class SingularSystemC : public NonlinearSystem
{
public:
  SingularSystemC() : NonlinearSystem( "Singular System C", SINGULAR_SYSTEM_BIBTEX, 2 ) {}

  virtual void evaluate( Vector const & X, Vector & f ) const override
  {
    real_type x = X( 0 );
    real_type y = X( 1 );
    f( 0 )      = exp( x ) - y - 1;
    f( 1 )      = x - y;
  }

  virtual void jacobian( Vector const & X, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type x      = X( 0 );
    J.insert( 0, 0 ) = exp( x );
    J.insert( 0, 1 ) = -1;
    J.insert( 1, 0 ) = 1;
    J.insert( 1, 1 ) = -1;
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
    x_vec.resize( 2 );
    auto & x0{ x_vec[0] };
    auto & x1{ x_vec[1] };
    x0.resize( n );
    x1.resize( n );
    x0 << 0.7, 0.7;
    x1 << 2, 1;
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class SingularSystemD : public NonlinearSystem
{
public:
  SingularSystemD() : NonlinearSystem( "Singular System D", SINGULAR_SYSTEM_BIBTEX, 2 ) {}

  virtual void evaluate( Vector const & X, Vector & f ) const override
  {
    real_type x = X( 0 );
    real_type y = X( 1 );
    f( 0 )      = power4( 6 * x - y );
    f( 1 )      = cos( x ) - 1 + y;
  }

  virtual void jacobian( Vector const & X, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type x      = X( 0 );
    real_type y      = X( 1 );
    J.insert( 0, 0 ) = 24 * power3( 6 * x - y );
    J.insert( 0, 1 ) = -4 * power3( 6 * x - y );
    J.insert( 1, 0 ) = -sin( x );
    J.insert( 1, 1 ) = 1;
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
    x_vec.resize( 2 );
    auto & x0{ x_vec[0] };
    auto & x1{ x_vec[1] };
    x0.resize( n );
    x1.resize( n );
    x0 << -0.5, -0.5;
    x1 << 4, 1;
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class SingularSystemE : public NonlinearSystem
{
public:
  SingularSystemE() : NonlinearSystem( "Singular System E", SINGULAR_SYSTEM_BIBTEX, 3 ) {}

  virtual void evaluate( Vector const & X, Vector & f ) const override
  {
    real_type x = X( 0 );
    real_type y = X( 1 );
    real_type z = X( 2 );
    f( 0 )      = 3 * x - cos( y * z ) - 0.5;
    f( 1 )      = x * x - 635 * y * y - 0.25;
    f( 2 )      = exp( -x * y ) + 20 * z + ( 10 * m_pi - 3 ) / 3;
  }

  virtual void jacobian( Vector const & X, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type x      = X( 0 );
    real_type y      = X( 1 );
    real_type z      = X( 2 );
    J.insert( 0, 0 ) = 3;
    J.insert( 0, 1 ) = z * sin( y * z );
    J.insert( 0, 2 ) = y * sin( y * z );
    J.insert( 1, 0 ) = 2 * x;
    J.insert( 1, 1 ) = -1270 * y;
    J.insert( 2, 0 ) = -y * exp( -x * y );
    J.insert( 2, 1 ) = -x * exp( -x * y );
    J.insert( 2, 2 ) = 20;
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 0.5, 0, -m_pi / 6;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 2 );
    auto & x0{ x_vec[0] };
    auto & x1{ x_vec[1] };
    x0.resize( n );
    x1.resize( n );

    x0 << 0.2, 0.2, -0.2;
    x1 << 1, 2, 1;
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class SingularSystemF : public NonlinearSystem
{
public:
  SingularSystemF() : NonlinearSystem( "Singular System F", SINGULAR_SYSTEM_BIBTEX, 3 ) {}

  virtual void evaluate( Vector const & X, Vector & f ) const override
  {
    real_type x = X( 0 );
    real_type y = X( 1 );
    real_type z = X( 2 );
    f( 0 )      = exp( x * x ) - 8 * x * sin( y );
    f( 1 )      = x + y - 1;
    f( 2 )      = power3( z - 1 );
  }

  virtual void jacobian( Vector const & X, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type x      = X( 0 );
    real_type y      = X( 1 );
    real_type z      = X( 2 );
    J.insert( 0, 0 ) = 2 * x * exp( x * x ) - 8 * sin( y );
    J.insert( 0, 1 ) = -8 * x * cos( y );
    J.insert( 0, 2 ) = 0;
    J.insert( 1, 0 ) = 1;
    J.insert( 1, 1 ) = 1;
    J.insert( 1, 2 ) = 0;
    J.insert( 2, 0 ) = 0;
    J.insert( 2, 1 ) = 0;
    J.insert( 2, 2 ) = 3 * power2( z - 1 );
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 0.175599, 0.824401, 1;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 2 );
    auto & x0{ x_vec[0] };
    auto & x1{ x_vec[1] };
    x0.resize( n );
    x1.resize( n );

    x0 << 0, 1, 2;
    x1 << 1, 1, 3;
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class SingularSystemP2 : public NonlinearSystem
{
public:
  SingularSystemP2() : NonlinearSystem( "Singular System Problem 2 (Ishihara, K. 2001)", SINGULAR_SYSTEM_BIBTEX, 3 ) {}

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = 4 * x( 0 ) - 2 * x( 1 ) + power2( x( 0 ) ) - 3;
    f( 1 ) = -x( 0 ) + 4 * x( 1 ) - x( 2 ) + power2( x( 1 ) ) - 3;
    f( 2 ) = -2 * x( 1 ) + 4 * x( 2 ) + power2( x( 2 ) ) - 3;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    J.insert( 0, 0 ) = 4 + 2 * x( 0 );
    J.insert( 0, 1 ) = -2;
    J.insert( 0, 2 ) = 0;

    J.insert( 1, 0 ) = -1;
    J.insert( 1, 1 ) = 4 + 2 * x( 1 );
    J.insert( 1, 2 ) = -1;

    J.insert( 2, 0 ) = 0;
    J.insert( 2, 1 ) = -2;
    J.insert( 2, 2 ) = 4 + 2 * x( 2 );

    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 1 );
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 5 );
    auto & x0{ x_vec[0] };
    auto & x1{ x_vec[1] };
    auto & x2{ x_vec[2] };
    auto & x3{ x_vec[3] };
    auto & x4{ x_vec[4] };
    x0.resize( n );
    x1.resize( n );
    x2.resize( n );
    x3.resize( n );
    x4.resize( n );

    x0 << -1.5, 0, -1.5;
    x1 << 4, 0, 4;
    x2 << -1, 5, -1;
    x3 << 4, 4, 4;
    x4 << -10, 0, 10;
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class SingularSystemP3 : public NonlinearSystem
{
public:
  SingularSystemP3() : NonlinearSystem( "Singular System Problem 3", SINGULAR_SYSTEM_BIBTEX, 2 ) {}

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = 2 / ( 1 + power2( x( 0 ) ) ) + sin( x( 1 ) - 1 ) - 1;
    f( 1 ) = 2 / ( 1 + power2( x( 1 ) ) ) + sin( x( 1 ) - 1 ) - 1;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    J.insert( 0, 0 ) = -4 * x( 0 ) / power2( power2( x( 0 ) ) + 1 );
    J.insert( 0, 1 ) = cos( x( 1 ) - 1 );
    J.insert( 1, 0 ) = 0;
    J.insert( 1, 1 ) = cos( x( 1 ) - 1 ) - 4 * x( 1 ) / power2( power2( x( 1 ) ) + 1 );
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 1 );
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 3 );
    auto & x0{ x_vec[0] };
    auto & x1{ x_vec[1] };
    auto & x2{ x_vec[2] };
    x0.resize( n );
    x1.resize( n );
    x2.resize( n );

    x0 << 0.5, 0.5;
    x1 << 2, 2;
    x2 << 0.1, 0.1;
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class SingularSystemP4 : public NonlinearSystem
{
public:
  SingularSystemP4() : NonlinearSystem( "Singular System Problem 4", SINGULAR_SYSTEM_BIBTEX, 2 ) {}

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = 1 + tan( 2 - 2 * cos( x( 0 ) ) ) - exp( sin( x( 0 ) ) );
    f( 1 ) = 1 + tan( 2 - 2 * cos( x( 1 ) ) ) - exp( sin( x( 1 ) ) );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    J.insert( 0, 0 ) = -cos( x( 0 ) ) * exp( sin( x( 0 ) ) ) +
                       2 * sin( x( 0 ) ) / power2( cos( -2 + 2 * cos( x( 0 ) ) ) );
    J.insert( 1, 1 ) = -cos( x( 1 ) ) * exp( sin( x( 1 ) ) ) +
                       2 * sin( x( 1 ) ) / power2( cos( -2 + 2 * cos( x( 1 ) ) ) );
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
    x_vec.resize( 3 );
    auto & x0{ x_vec[0] };
    auto & x1{ x_vec[1] };
    auto & x2{ x_vec[2] };
    x0.resize( n );
    x1.resize( n );
    x2.resize( n );
    x0 << 3, 0;
    x1 << 0, 0.5;
    x2 << -0.5, -0.5;
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class SingularSystemP5 : public NonlinearSystem
{
public:
  SingularSystemP5() : NonlinearSystem( "Singular System Problem 5", SINGULAR_SYSTEM_BIBTEX, 2 ) {}

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = exp( x( 0 ) ) + x( 1 ) - 1;
    f( 1 ) = exp( x( 1 ) ) + x( 0 ) - 1;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    J.insert( 0, 0 ) = exp( x( 0 ) );
    J.insert( 0, 1 ) = 1;
    J.insert( 1, 0 ) = 1;
    J.insert( 1, 1 ) = exp( x( 1 ) );

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
    x0 << -0.5, -0.5;
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class SingularSystemP6 : public NonlinearSystem
{
public:
  SingularSystemP6() : NonlinearSystem( "Singular System Problem 6", SINGULAR_SYSTEM_BIBTEX, 3 ) {}

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = cos( x( 0 ) ) - 9 + 3 * x( 0 ) + 8 * exp( x( 1 ) );
    f( 1 ) = cos( x( 1 ) ) - 9 + 3 * x( 1 ) + 8 * exp( x( 0 ) );
    f( 2 ) = cos( x( 2 ) ) - x( 2 ) - 1;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    J.insert( 0, 0 ) = -sin( x( 0 ) ) + 3;
    J.insert( 0, 1 ) = 8 * exp( x( 1 ) );
    J.insert( 1, 0 ) = 8 * exp( x( 0 ) );
    J.insert( 1, 1 ) = -sin( x( 1 ) ) + 3;
    J.insert( 2, 2 ) = -sin( x( 2 ) ) - 1;
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
    x_vec.resize( 4 );
    auto & x0{ x_vec[0] };
    auto & x1{ x_vec[1] };
    auto & x2{ x_vec[2] };
    auto & x3{ x_vec[3] };
    x0.resize( n );
    x1.resize( n );
    x2.resize( n );
    x3.resize( n );

    x0 << -1, -1, -1;
    x1 << 3, 3, 3;
    x2 << 0.5, 0.5, 0.5;
    x3 << -3, -3, -3;
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class SingularSystemP7 : public NonlinearSystem
{
public:
  SingularSystemP7() : NonlinearSystem( "Singular System Problem 7 (Ishihara, K. 2001)", SINGULAR_SYSTEM_BIBTEX, 2 ) {}

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = 4 * x( 0 ) - 2 * x( 1 ) + power2( x( 0 ) ) - 3;
    f( 1 ) = -2 * x( 0 ) + 4 * x( 1 ) + power2( x( 0 ) ) - 3;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    J.insert( 0, 0 ) = 4 + 2 * x( 0 );
    J.insert( 0, 1 ) = -2;
    J.insert( 1, 0 ) = -2 + 2 * x( 0 );
    J.insert( 1, 1 ) = 4;
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 1 );
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 4 );
    auto & x0{ x_vec[0] };
    auto & x1{ x_vec[1] };
    auto & x2{ x_vec[2] };
    auto & x3{ x_vec[3] };
    x0.resize( n );
    x1.resize( n );
    x2.resize( n );
    x3.resize( n );

    x0 << 3, 3;
    x1 << 0, -1.5;
    x2 << -2, 3;
    x3 << 0, 2;
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class SingularSystemP8 : public NonlinearSystem
{
public:
  SingularSystemP8() : NonlinearSystem( "Singular System Problem 8", SINGULAR_SYSTEM_BIBTEX, 2 ) {}

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = sqrt( 3.0 ) * power2( x( 0 ) ) - power2( x( 1 ) );
    f( 1 ) = cos( x( 0 ) ) - 1 / ( 1 + power2( x( 1 ) ) );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    J.insert( 0, 0 ) = 2 * sqrt( 3.0 ) * x( 0 );
    J.insert( 0, 1 ) = -2 * x( 1 );
    J.insert( 1, 0 ) = -sin( x( 0 ) );
    J.insert( 1, 1 ) = 2 * x( 1 ) / power2( power2( x( 1 ) ) + 1 );
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
    x0 << 0.5, 1;
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class SingularSystemP9 : public NonlinearSystem
{
public:
  SingularSystemP9() : NonlinearSystem( "Singular System Problem 9", SINGULAR_SYSTEM_BIBTEX, 2 ) {}

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = power2( x( 0 ) ) - power2( x( 1 ) );
    f( 1 ) = 3 * power2( x( 0 ) ) - 3 * power2( x( 1 ) );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    J.insert( 0, 0 ) = 2 * x( 0 );
    J.insert( 0, 1 ) = -2 * x( 1 );
    J.insert( 1, 0 ) = 6 * x( 0 );
    J.insert( 1, 1 ) = -6 * x( 1 );
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
    x_vec.resize( 4 );
    auto & x0{ x_vec[0] };
    auto & x1{ x_vec[1] };
    auto & x2{ x_vec[2] };
    auto & x3{ x_vec[3] };
    x0.resize( n );
    x1.resize( n );
    x2.resize( n );
    x3.resize( n );

    x0 << 0.5, 0.4;
    x1 << -0.5, -0.4;
    x2 << 0.3, -0.5;
    x3 << 0.4, 0.5;
  }
};
