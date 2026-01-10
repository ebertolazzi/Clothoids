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

class BroydenBandedFunction : public NonlinearSystem
{
  integer const ml;
  integer const mu;

public:
  BroydenBandedFunction()
    : NonlinearSystem(
        "Broyden Banded Function",
        "@article{Broyden:1971,\n"
        "  author  = {Broyden, C. G.},\n"
        "  title   = {The convergence of an algorithm for solving sparse "
        "nonlinear systems},\n"
        "  journal = {Mathematics of Computation},\n"
        "  volume  = {25},\n"
        "  year    = {1971},\n"
        "  pages   = {285--294},\n"
        "  doi     = {10.2307/2004922},\n"
        "}\n\n"
        "@article{More:1981,\n"
        "  author  = {Mor{\'e}, Jorge J. and Garbow, Burton S. and "
        "Hillstrom, Kenneth E.},\n"
        "  title   = {Testing Unconstrained Optimization Software},\n"
        "  journal = {ACM Trans. Math. Softw.},\n"
        "  volume  = {7},\n"
        "  number  = {1},\n"
        "  year    = {1981},\n"
        "  pages   = {17--41},\n"
        "  doi     = {10.1145/355934.355936},\n"
        "}\n",
        10 )
    , ml( 5 )
    , mu( 1 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    for ( integer k = 0; k < n; ++k )
    {
      integer   k1   = max( static_cast<integer>( 0 ), k - ml );
      integer   k2   = min( n - 1, k + mu );
      real_type temp = 0;
      for ( integer j = k1; j <= k2; ++j )
      {
        if ( j != k ) temp += x( j ) * ( 1 + x( j ) );
      }
      f( k ) = x( k ) * ( 2 + 5 * x( k ) * x( k ) ) + 1 - temp;
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer k = 0; k < n; ++k )
    {
      integer k1 = max( static_cast<integer>( 0 ), k - ml );
      integer k2 = min( n - 1, k + mu );
      for ( integer j = k1; j <= k2; ++j )
      {
        if ( j != k ) J.insert( k, j ) = -( 1 + 2 * x( j ) );
      }
      J.insert( k, k ) = 2 + 15 * x( k ) * x( k );
    }
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
