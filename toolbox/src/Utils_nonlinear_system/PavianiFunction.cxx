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

class PavianiFunction : public NonlinearSystem
{
public:
  // sum log(xi-2)^2+log(xi-10)^2 - prod( xi) ^(1/5)
  PavianiFunction()
    : NonlinearSystem(
        "Paviani function",
        "@book{himmelblau:1972,\n"
        "  author    = {Himmelblau, D.M.},\n"
        "  title     = {Applied nonlinear programming},\n"
        "  year      = {1972},\n"
        "  publisher = {McGraw-Hill}\n"
        "}\n",
        10 )
  {
  }

  real_type loglog( real_type x ) const { return power2( log( 10 - x ) ) + power2( log( x - 2 ) ); }

  real_type loglog_D( real_type x ) const
  {
    real_type t1 = 10 - x;
    real_type t2 = x - 2;
    return 2 * ( log( t2 ) / t2 - log( t1 ) / t1 );
  }

  real_type loglog_DD( real_type x ) const
  {
    real_type t1 = 10 - x;
    real_type t2 = x - 2;
    return 2 * ( ( 1 - log( t1 ) ) / ( t1 * t1 ) + ( 1 - log( t2 ) ) / ( t2 * t2 ) );
  }

  real_type mul( Vector const & x ) const
  {
    real_type res = 1;
    for ( integer i = 0; i < 10; ++i ) res *= power2( std::abs( x( i ) ) );
    return res;
  }

  real_type mul_D( Vector const & x, integer k ) const
  {
    real_type res = 1;
    for ( integer i = 0; i < 10; ++i )
    {
      res *= x( i );
      if ( i == k )
        res *= 2;
      else
        res *= x( i );
    }
    return res;
  }

  real_type mul_DD( Vector const & x, integer i, integer j ) const
  {
    real_type res = 1;
    if ( i == j )
    {
      for ( integer k = 0; k < 10; ++k )
      {
        if ( i == k )
          res *= 2;
        else
          res *= power2( x( k ) );
      }
    }
    else
    {
      for ( integer k = 0; k < 10; ++k )
      {
        res *= x( k );
        if ( i == k || j == k )
          res *= 2;
        else
          res *= x( k );
      }
    }
    return res;
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    for ( integer i = 0; i < 10; ++i )
    {
      if ( x( i ) > 2 && x( i ) < 10 )
        f( i ) = loglog_D( x( i ) ) - mul_D( x, i );
      else
        f( i ) = nan( "PavianiFunction" );
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer i = 0; i < 10; ++i )
    {
      for ( integer j = 0; j < 10; ++j )
      {
        real_type tmp = -mul_DD( x, i, j );
        if ( i == j ) tmp += loglog_DD( x( i ) );
        J.insert( i, j ) = tmp;
      }
    }
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    /*x(0) = 2.001;
    x(1) = 9.999;
    x(2) = 2.001;
    x(3) = 9.999;
    x(4) = 2.001;
    x(5) = 9.999;
    x(6) = 2.001;
    x(7) = 9.999;
    x(8) = 2.001;
    x(9) = 9.999;
    */
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 2.1, 9.9, 2.1, 9.9, 2.1, 9.9, 2.1, 9.9, 2.1, 9.9;
  }

  virtual void check_if_admissible( Vector const & x ) const override
  {
    for ( integer i = 0; i < 10; ++i )
    {
      UTILS_ASSERT( x( i ) > 2 && x( i ) < 10, "x[{}] = {} must be in (2,10)", i, x( i ) );
    }
  }

  virtual void bounding_box( Vector & L, Vector & U ) const override
  {
    U.fill( 10 );
    L.fill( 2 );
  }
};
