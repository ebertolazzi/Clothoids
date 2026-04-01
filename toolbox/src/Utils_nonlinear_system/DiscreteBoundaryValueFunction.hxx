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

class DiscreteBoundaryValueFunction : public NonlinearSystem
{
  real_type h;

public:
  DiscreteBoundaryValueFunction( integer neq )
    : NonlinearSystem(
        "Discrete boundary value function",
        "@article{More:1979,\n"
        "  author  = {Mor{\'e}, Jorge J. and Cosnard, Michel Y.},\n"
        "  title   = {Numerical Solution of Nonlinear Equations},\n"
        "  journal = {ACM Trans. Math. Softw.},\n"
        "  year    = {1979},\n"
        "  volume  = {5},\n"
        "  number  = {1},\n"
        "  pages   = {64--85},\n"
        "  doi     = {10.1145/355815.355820},\n"
        "}\n\n"
        "@article{More:1981,\n"
        "  author  = {Mor{\'e}, Jorge J. and Garbow, Burton S. and "
        "Hillstrom, Kenneth E.},\n"
        "  title   = {Testing Unconstrained Optimization Software},\n"
        "  journal = {ACM Trans. Math. Softw.},\n"
        "  year    = {1981},\n"
        "  volume  = {7},\n"
        "  number  = {1},\n"
        "  pages   = {17--41},\n"
        "  doi     = {10.1145/355934.355936},\n"
        "}\n",
        neq )
    , h( 1 / real_type( n + 1 ) )
  {
    check_min_equations( n, 1 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    for ( integer k = 0; k < n; ++k )
    {
      f( k ) = 2 * x( k ) + 0.5 * power2( h ) * power3( ( x( k ) + 1 ) + ( k + 1 ) * h );
      if ( k > 0 ) f( k ) -= x( k - 1 );
      if ( k < n - 1 ) f( k ) -= x( k + 1 );
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer i = 0; i < n; ++i ) J.insert( i, i ) = 2 + 1.5 * h * h * power2( x( i ) + h * ( i + 1 ) + 1 );

    for ( integer i = 0; i < n - 1; ++i )
    {
      J.insert( i + 1, i ) = -1;
      J.insert( i, i + 1 ) = -1;
    }
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.clear();
    if ( n == 2 )
    {
      x_vec.resize( 1 );
      auto & x0{ x_vec[0] };
      x0.resize( n );
      x0 << -0.128246763033732, -0.159267567244641;
    }
    else if ( n == 5 )
    {
      x_vec.resize( 1 );
      auto & x0{ x_vec[0] };
      x0.resize( n );
      x0 << -0.0750221292923205, -0.131976210352191, -0.164848771909337, -0.164664680215801, -0.117417651684194;
    }
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    for ( integer k = 0; k < n; ++k ) x0( k ) = real_type( ( k + 1 ) * ( k - n ) ) / power2( n + 1.0 );
  }
};
