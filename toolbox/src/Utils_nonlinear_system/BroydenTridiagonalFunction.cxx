/*\
 |
 |  Author:
 |    Enrico Bertolazzi
 |    University of Trento
 |    Department of Industrial Engineering
 |    Via Sommarive 9, I-38123, Povo, Trento, Italy
 |    email: enrico.bertolazzi@unitn.it
\*/

// TEST 222

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class BroydenTridiagonalFunction : public NonlinearSystem
{
  real_type const alpha;
  real_type const beta;

public:
  BroydenTridiagonalFunction( real_type _alpha, real_type _beta, integer _neq )
    : NonlinearSystem(
        "Broyden tridiagonal function",
        "@article{Broyden:1965,\n"
        "  author  = {Broyden, C. G.},\n"
        "  title   = {A class of methods for solving nonlinear "
        "simultaneous equations},\n"
        "  journal = {Mathematics of Computation},\n"
        "  volume  = {19},\n"
        "  year    = {1965},\n"
        "  pages   = {577--593},\n"
        "  doi     = {10.2307/2003941}\n"
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
        _neq )
    , alpha( _alpha )
    , beta( _beta )
  {
    check_min_equations( n, 1 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    for ( integer k = 0; k < n; ++k )
    {
      f( k ) = ( 3 - alpha * x( k ) ) * x( k ) + beta;
      if ( k > 0 ) f( k ) -= x( k - 1 );
      if ( k < n - 1 ) f( k ) -= 2 * x( k + 1 );
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer k = 0; k < n; ++k ) J.insert( k, k ) = 3 - 2 * alpha * x( k );
    for ( integer k = 1; k < n; ++k ) J.insert( k, k - 1 ) = -1;
    for ( integer k = 0; k < n - 1; ++k ) J.insert( k, k + 1 ) = -2;
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( -1 );
  }
};
