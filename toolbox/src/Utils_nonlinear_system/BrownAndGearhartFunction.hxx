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

class BrownAndGearhartFunction : public NonlinearSystem
{
  real_type const sqrt2;

public:
  BrownAndGearhartFunction()
    : NonlinearSystem(
        "Brown and Gearhart function",
        "@Article{Brow1971,\n"
        "  author  = {Brow, Kenneth M. and Gearhart, William B.},\n"
        "  title   = {Deflation techniques for the calculation of further "
        "solutions\n"
        "             of a nonlinear system},\n"
        "  journal = {Numerische Mathematik},\n"
        "  year    = {1971},\n"
        "  volume  = {16},\n"
        "  number  = {4},\n"
        "  pages   = {334--342},\n"
        "  doi     = \"10.1007/BF02165004\",\n"
        "}\n",
        3 )
    , sqrt2( sqrt( 2.0 ) )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = power2( x( 0 ) ) + 2 * power2( x( 1 ) ) - 4;
    f( 1 ) = power2( x( 0 ) ) + power2( x( 1 ) ) + x( 2 ) - 8;
    f( 2 ) = power2( x( 0 ) - 1 ) + power2( 2 * x( 1 ) - sqrt2 ) + power2( x( 2 ) - 5 ) - 4;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    J.insert( 0, 0 ) = 2 * x( 0 );
    J.insert( 0, 1 ) = 4 * x( 1 );
    J.insert( 0, 2 ) = 0;

    J.insert( 1, 0 ) = 2 * x( 0 );
    J.insert( 1, 1 ) = 2 * x( 1 );
    J.insert( 1, 2 ) = 1;

    J.insert( 2, 0 ) = 2 * ( x( 0 ) - 1 );
    J.insert( 2, 1 ) = 4 * ( 2 * x( 1 ) - sqrt2 );
    J.insert( 2, 2 ) = 2 * ( x( 2 ) - 5 );

    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 0, sqrt2, 6;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 1, 0.7, 5;
  }
};
