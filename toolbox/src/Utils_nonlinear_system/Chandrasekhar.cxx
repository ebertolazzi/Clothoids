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

class Chandrasekhar : public NonlinearSystem
{
  using Matrix = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;

  Vector          mu;
  real_type const w;

public:
  Chandrasekhar( real_type c, integer neq )
    : NonlinearSystem(
        "Chandrasekhar function",
        "@book{Kelley:1995,\n"
        "  author    = {Kelley, C.},\n"
        "  title     = {Iterative Methods for Linear and Nonlinear "
        "Equations},\n"
        "  publisher = {Society for Industrial and Applied Mathematics},\n"
        "  year      = {1995},\n"
        "  doi       = {10.1137/1.9781611970944},\n"
        "}\n\n"
        "@book{chandrasekhar1960,\n"
        "  author    = {Chandrasekhar, S.},\n"
        "  title     = {Radiative Transfer},\n"
        "  year      = {1960},\n"
        "  series    = {Dover Books on Intermediate and Advanced "
        "Mathematics},\n"
        "  publisher = {Dover Publications},\n"
        "  isbn      = {9780486605906}\n"
        "}\n",
        neq )
    , w( c / ( 2 * neq ) )
  {
    mu.resize( neq );
    for ( integer i = 0; i < neq; ++i ) mu( i ) = i + 0.5;
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    for ( integer i = 0; i < n; ++i )
    {
      real_type tmp = 0;
      for ( integer j = 0; j < n; ++j ) tmp += mu( j ) * x( j ) / ( mu( i ) + mu( j ) );
      f( i ) = x( i ) - 1 / ( 1 - w * tmp );
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    Matrix J_full( n, n );
    J_full.setZero();
    for ( integer i = 0; i < n; ++i )
    {
      real_type tmp = 0;
      for ( integer j = 0; j < n; ++j ) tmp += mu( j ) * x( j ) / ( mu( i ) + mu( j ) );
      tmp = -w / power2( 1 - w * tmp );
      for ( integer j = 0; j < n; ++j )
      {
        J_full( i, j ) = tmp * mu( j ) / ( mu( i ) + mu( j ) );
        if ( i == j ) J_full( i, j ) += 1;
      }
    }
    J.resize( n, n );
    J = J_full.sparseView();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 10 );
  }
};
