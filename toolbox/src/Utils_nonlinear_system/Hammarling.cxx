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

#define HAMMARLING_BIBTEX                          \
  "n by n matrix square root problem (hammarling)" \
  "S. J. Hammarling, private communication to P. E. Gill.\n"

class Hammarling2x2matrixSquareRoot : public NonlinearSystem
{
  real_type const a00;
  real_type const a01;
  real_type const a10;
  real_type const a11;

public:
  Hammarling2x2matrixSquareRoot()
    : NonlinearSystem( "Hammarling 2 by 2 matrix square root problem", HAMMARLING_BIBTEX, 4 )
    , a00( 1.0E-4 )
    , a01( 1 )
    , a10( 0 )
    , a11( 1E-4 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = ( x( 0 ) * x( 0 ) + x( 1 ) * x( 2 ) ) - a00;
    f( 1 ) = ( x( 0 ) * x( 1 ) + x( 1 ) * x( 3 ) ) - a01;
    f( 2 ) = ( x( 2 ) * x( 0 ) + x( 3 ) * x( 2 ) ) - a10;
    f( 3 ) = ( x( 2 ) * x( 1 ) + x( 3 ) * x( 3 ) ) - a11;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    J.insert( 0, 0 ) = 2 * x( 0 );
    J.insert( 0, 1 ) = x( 2 );
    J.insert( 0, 2 ) = x( 1 );
    J.insert( 0, 3 ) = 0;

    J.insert( 1, 0 ) = x( 1 );
    J.insert( 1, 1 ) = x( 0 ) + x( 3 );
    J.insert( 1, 2 ) = 0;
    J.insert( 1, 3 ) = x( 1 );

    J.insert( 2, 0 ) = x( 2 );
    J.insert( 2, 1 ) = 0;
    J.insert( 2, 2 ) = x( 0 ) + x( 3 );
    J.insert( 2, 3 ) = x( 2 );

    J.insert( 3, 0 ) = 0;
    J.insert( 3, 1 ) = x( 2 );
    J.insert( 3, 2 ) = x( 1 );
    J.insert( 3, 3 ) = 2 * x( 3 );

    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 1E-2, 5E+1, 0, 1E-2;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 1, 0, 0, 1;
  }
};

/*
 * Hammarling 3 by 3 matrix square root problem
 */

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class Hammarling3x3matrixSquareRootProblem : public NonlinearSystem
{
  real_type a00, a01, a02, a10, a11, a12, a20, a21, a22;

  Vector xe;

  real_type f0( Vector const & x ) const { return x( 0 ) * x( 0 ) + x( 1 ) * x( 3 ) + x( 2 ) * x( 6 ); }

  real_type f1( Vector const & x ) const { return x( 0 ) * x( 1 ) + x( 1 ) * x( 4 ) + x( 2 ) * x( 7 ); }

  real_type f2( Vector const & x ) const { return x( 0 ) * x( 2 ) + x( 1 ) * x( 5 ) + x( 2 ) * x( 8 ); }

  real_type f3( Vector const & x ) const { return x( 3 ) * x( 0 ) + x( 4 ) * x( 3 ) + x( 5 ) * x( 6 ); }

  real_type f4( Vector const & x ) const { return x( 3 ) * x( 1 ) + x( 4 ) * x( 4 ) + x( 5 ) * x( 7 ); }

  real_type f5( Vector const & x ) const { return x( 3 ) * x( 2 ) + x( 4 ) * x( 5 ) + x( 5 ) * x( 8 ); }

  real_type f6( Vector const & x ) const { return x( 6 ) * x( 0 ) + x( 7 ) * x( 3 ) + x( 8 ) * x( 6 ); }

  real_type f7( Vector const & x ) const { return x( 6 ) * x( 1 ) + x( 7 ) * x( 4 ) + x( 8 ) * x( 7 ); }

  real_type f8( Vector const & x ) const { return x( 6 ) * x( 2 ) + x( 7 ) * x( 5 ) + x( 8 ) * x( 8 ); }

public:
  Hammarling3x3matrixSquareRootProblem(
    string const & n,
    real_type      x0,
    real_type      x1,
    real_type      x2,
    real_type      x3,
    real_type      x4,
    real_type      x5,
    real_type      x6,
    real_type      x7,
    real_type      x8 )
    : NonlinearSystem( n, HAMMARLING_BIBTEX, 9 )
  {
    xe.resize( 9 );

    xe[0] = x0;
    xe[1] = x1;
    xe[2] = x2;
    xe[3] = x3;
    xe[4] = x4;
    xe[5] = x5;
    xe[6] = x6;
    xe[7] = x7;
    xe[8] = x8;

    a00 = f0( xe );
    a01 = f1( xe );
    a02 = f2( xe );
    a10 = f3( xe );
    a11 = f4( xe );
    a12 = f5( xe );
    a20 = f6( xe );
    a21 = f7( xe );
    a22 = f8( xe );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = f0( x ) - a00;
    f( 1 ) = f1( x ) - a01;
    f( 2 ) = f2( x ) - a02;
    f( 3 ) = f3( x ) - a10;
    f( 4 ) = f4( x ) - a11;
    f( 5 ) = f5( x ) - a12;
    f( 6 ) = f6( x ) - a20;
    f( 7 ) = f7( x ) - a21;
    f( 8 ) = f8( x ) - a22;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    J.insert( 0, 0 ) = 2 * x( 0 );
    J.insert( 0, 1 ) = x( 3 );
    J.insert( 0, 2 ) = x( 6 );
    J.insert( 0, 3 ) = x( 1 );
    J.insert( 0, 6 ) = x( 2 );

    J.insert( 1, 0 ) = x( 1 );
    J.insert( 1, 1 ) = x( 0 ) + x( 4 );
    J.insert( 1, 2 ) = x( 7 );
    J.insert( 1, 4 ) = x( 1 );
    J.insert( 1, 7 ) = x( 2 );

    J.insert( 2, 0 ) = x( 2 );
    J.insert( 2, 1 ) = x( 5 );
    J.insert( 2, 2 ) = x( 0 ) + x( 8 );
    J.insert( 2, 5 ) = x( 1 );
    J.insert( 2, 8 ) = x( 2 );

    J.insert( 3, 0 ) = x( 3 );
    J.insert( 3, 3 ) = x( 0 ) + x( 4 );
    J.insert( 3, 4 ) = x( 3 );
    J.insert( 3, 5 ) = x( 6 );
    J.insert( 3, 6 ) = x( 5 );

    J.insert( 4, 1 ) = x( 3 );
    J.insert( 4, 3 ) = x( 1 );
    J.insert( 4, 4 ) = 2 * x( 4 );
    J.insert( 4, 5 ) = x( 7 );
    J.insert( 4, 7 ) = x( 5 );

    J.insert( 5, 2 ) = x( 3 );
    J.insert( 5, 3 ) = x( 2 );
    J.insert( 5, 4 ) = x( 5 );
    J.insert( 5, 5 ) = x( 4 ) + x( 8 );
    J.insert( 5, 8 ) = x( 5 );

    J.insert( 6, 0 ) = x( 6 );
    J.insert( 6, 3 ) = x( 7 );
    J.insert( 6, 6 ) = x( 0 ) + x( 8 );
    J.insert( 6, 7 ) = x( 3 );
    J.insert( 6, 8 ) = x( 6 );

    J.insert( 7, 1 ) = x( 6 );
    J.insert( 7, 4 ) = x( 7 );
    J.insert( 7, 6 ) = x( 1 );
    J.insert( 7, 7 ) = x( 4 ) + x( 8 );
    J.insert( 7, 8 ) = x( 7 );

    J.insert( 8, 2 ) = x( 6 );
    J.insert( 8, 5 ) = x( 7 );
    J.insert( 8, 6 ) = x( 2 );
    J.insert( 8, 7 ) = x( 5 );
    J.insert( 8, 8 ) = 2 * x( 8 );

    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 = xe;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 1, 0, 0, 0, 1, 0, 0, 0, 1;
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class Hammarling3x3matrixSquareRootProblemN1 : public Hammarling3x3matrixSquareRootProblem
{
public:
  Hammarling3x3matrixSquareRootProblemN1()
    : Hammarling3x3matrixSquareRootProblem(
        "Hammarling 3 by 3 matrix square root problem N.1",
        0.01,
        50,
        0,
        0,
        0.01,
        0,
        0,
        0,
        0.01 )
  {
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class Hammarling3x3matrixSquareRootProblemN2 : public Hammarling3x3matrixSquareRootProblem
{
public:
  Hammarling3x3matrixSquareRootProblemN2()
    : Hammarling3x3matrixSquareRootProblem(
        "Hammarling 3 by 3 matrix square root problem N.2",
        0,
        0,
        1,
        1,
        1,
        0,
        0,
        1,
        0 )
  {
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class Hammarling3x3matrixSquareRootProblemN3 : public Hammarling3x3matrixSquareRootProblem
{
public:
  Hammarling3x3matrixSquareRootProblemN3()
    : Hammarling3x3matrixSquareRootProblem(
        "Hammarling 3 by 3 matrix square root problem N.3",
        1,
        1,
        1,
        0,
        0,
        0,
        0,
        0,
        0 )
  {
  }
};
