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

class GriewankFunction : public NonlinearSystem
{
public:
  GriewankFunction( integer n )
    : NonlinearSystem(
        "Griewank function",
        "@Article{Griewan:k1981,\n"
        "  author  = {Griewank, A. O.},\n"
        "  title   = {Generalized descent for global optimization},\n"
        "  journal = {Journal of Optimization Theory and Applications},\n"
        "  year    = {1981},\n"
        "  volume  = {34},\n"
        "  number  = {1},\n"
        "  pages   = {11--39},\n"
        "  doi     = {10.1007/BF00933356}\n"
        "}\n",
        n )
  {
    UTILS_ASSERT( n >= 2 && n <= 20, "GriewankFunction(n={}) must be in range [2..20]", n );
  }

  real_type grad( Vector const & x, integer k ) const
  {
    real_type f = 1 / sqrt( k + 1.0 );
    for ( integer i = 0; i < n; ++i )
    {
      real_type t = x( i ) / sqrt( i + 1.0 );
      if ( i == k )
        f *= sin( t );
      else
        f *= cos( t );
    }
    return f;
  }

  real_type hess( Vector const & x, integer i, integer j ) const
  {
    if ( i == j )
    {
      real_type f = 1 / ( i + 1.0 );
      for ( integer k = 0; k < n; ++k )
      {
        real_type t = x( k ) / sqrt( k + 1.0 );
        f *= cos( t );
      }
      return f;
    }
    else
    {
      real_type f = -1 / sqrt( ( i + 1.0 ) * ( j + 1.0 ) );
      for ( integer k = 0; k < n; ++k )
      {
        real_type t = x( k ) / sqrt( k + 1.0 );
        if ( k == i || k == j )
          f *= sin( t );
        else
          f *= cos( t );
      }
      return f;
    }
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    for ( integer k = 0; k < n; ++k ) f( k ) = grad( x, k ) + x( k ) / 2000.0;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer i = 0; i < n; ++i )
    {
      for ( integer j = 0; j < n; ++j )
      {
        real_type tmp = hess( x, i, j );
        if ( i == j ) tmp += 1.0 / 2000.0;
        J.insert( i, j ) = tmp;
      }
    }
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 2 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    for ( integer i{ 0 }; i < n; i += 2 )
    {
      x0( i + 0 ) = -50.0;
      if ( i + 1 < n ) x0( i + 1 ) = 70.0;
    }
    auto & x1{ x_vec[1] };
    x1.resize( n );
    x1 = 10 * x0;
  }

  virtual void check_if_admissible( Vector const & x ) const override
  {
    for ( integer i = 0; i < n; ++i ) UTILS_ASSERT( std::abs( x( i ) ) < 1000, "Bad range" );
  }

  virtual void bounding_box( Vector & L, Vector & U ) const override
  {
    L.fill( -1000 );
    U.fill( 1000 );
  }
};
