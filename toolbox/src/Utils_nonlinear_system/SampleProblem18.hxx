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

class SampleProblem18 : public NonlinearSystem
{
  real_type const epsi;

public:
  SampleProblem18() : NonlinearSystem( "Sample problem 18", "no doc", 2 ), epsi( 1e-8 ) {}

  real_type fun( real_type x ) const
  {
    real_type x2 = x * x;
    if ( std::abs( x2 ) > epsi )
      return ( 1 - exp( -x2 ) ) / x;
    else
      return x * ( 1 + x2 * ( x2 / 6.0 - 0.5 ) );
  }

  real_type fun_1( real_type x ) const
  {
    real_type x2    = x * x;
    real_type expx2 = exp( -x2 );
    if ( std::abs( x2 ) > epsi )
      return 2 * expx2 - ( 1 - expx2 ) / x2;
    else
      return 1 + x2 * ( ( 5.0 / 6.0 ) * x2 - 1.5 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = power2( x( 1 ) ) * fun( x( 0 ) );
    f( 1 ) = x( 0 ) * fun( x( 1 ) );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    J.insert( 0, 0 ) = power2( x( 1 ) ) * fun_1( x( 0 ) );
    J.insert( 0, 1 ) = 2 * x( 1 ) * fun( x( 0 ) );
    J.insert( 1, 0 ) = fun( x( 1 ) );
    J.insert( 1, 1 ) = x( 0 ) * fun_1( x( 1 ) );
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
    x0 << 2, 2;
  }

  virtual void check_if_admissible( Vector const & x ) const override
  {
    for ( integer i{ 0 }; i < n; ++i ) UTILS_ASSERT( std::abs( x( i ) ) < 4, "Bad range" );
  }

  virtual void bounding_box( Vector & L, Vector & U ) const override
  {
    U.fill( 4 );
    L.fill( -4 );
  }
};
