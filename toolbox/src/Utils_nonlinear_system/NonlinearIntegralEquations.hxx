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

class NonlinearIntegralEquations : public NonlinearSystem
{
  using Matrix = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;

public:
  NonlinearIntegralEquations()
    : NonlinearSystem(
        "Nonlinear Integral Equations",
        "@article{More:1979,\n"
        "  author  = {Mor{\'e}, Jorge J. and Cosnard, Michel Y.},\n"
        "  title   = {Numerical Solution of Nonlinear Equations},\n"
        "  journal = {ACM Trans. Math. Softw.},\n"
        "  year    = {1979},\n"
        "  volume  = {5},\n"
        "  number  = {1},\n"
        "  pages   = {64--85},\n"
        "  doi     = {10.1145/355815.355820},\n"
        "}\n",
        100 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type h = 1.0 / ( n - 1.0 );

    for ( integer i = 0; i < n; ++i ) f( i ) = x( i );

    real_type si = 0;
    for ( integer i = 1; i < n - 1; ++i )
    {
      real_type t = h * i;
      si += power3( x( i ) + t + 1 );
      f( i ) += 0.5 * ( 1 - t ) * si;
    }

    si = 0;
    for ( integer i = n - 2; i > 0; --i )
    {
      real_type t = h * i;
      f( i ) += 0.5 * ( 1 - t ) * t * si;
      si += power2( x( i ) + t + 1 );
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    Matrix J_full( n, n );
    J_full.setZero();

    real_type h = 1.0 / ( n - 1.0 );

    for ( integer i = 0; i < n; ++i ) J_full( i, i ) = 1;

    for ( integer i = 0; i < n; ++i )
    {
      real_type t = h * i;
      for ( integer j = 1; j <= i; ++j ) { J_full( i, j ) += 1.5 * ( 1 - t ) * power2( x( j ) + h * j + 1 ); }
      for ( integer j = i + 1; j < n - 1; ++j ) { J_full( i, j ) += ( 1 - t ) * t * ( x( j ) + h * j + 1 ); }
    }
    J.resize( n, n );
    J = J_full.sparseView();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    real_type h = 1.0 / ( n - 1.0 );
    for ( integer i{ 0 }; i < n; ++i )
    {
      real_type t = h * i;
      x0( i )     = t * ( t - 1 );
    }
  }
};
