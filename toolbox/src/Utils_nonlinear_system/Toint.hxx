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

class Toint225 : public NonlinearSystem
{
public:
  Toint225( integer neq )
    : NonlinearSystem(
        "Toint N.225",
        "@Article{Spedicato1997,\n"
        "  author  = {Spedicato, E. and Huang, Z.},\n"
        "  title   = {Numerical experience with newton-like methods\n"
        "             for nonlinear algebraic systems},\n"
        "  journal = {Computing},\n"
        "  year    = {1997},\n"
        "  volume  = {58},\n"
        "  number  = {1},\n"
        "  pages   = {69--89},\n"
        "  doi     = {10.1007/BF02684472},\n"
        "}\n",
        neq )
  {
    check_min_equations( n, 2 );
  }

  real_type phi1( real_type s, real_type t ) const { return 3 * s * s + 2 * t - 5 + sin( s - t ) * sin( s + t ); }

  real_type phi2( real_type s, real_type t ) const { return 4 * t - 3 + s * exp( s - t ); }

  real_type phi1_1( real_type s, real_type t ) const
  {
    return 6 * s + cos( s - t ) * sin( s + t ) + sin( s - t ) * cos( s + t );
  }

  real_type phi1_2( real_type s, real_type t ) const
  {
    return 2 - cos( s - t ) * sin( s + t ) + sin( s - t ) * cos( s + t );
  }

  real_type phi2_1( real_type s, real_type t ) const { return ( 1 + s ) * exp( s - t ); }

  real_type phi2_2( real_type s, real_type t ) const { return 4 - s * exp( s - t ); }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = phi1( x( 0 ), x( 1 ) );
    for ( integer k = 1; k < n - 1; ++k ) f( k ) = phi1( x( k ), x( k + 1 ) ) + phi2( x( k - 1 ), x( k ) );
    f( n - 1 ) = phi2( x( n - 2 ), x( n - 1 ) );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    J.insert( 0, 0 ) = phi1_1( x( 0 ), x( 1 ) );
    J.insert( 0, 1 ) = phi1_2( x( 0 ), x( 1 ) );
    for ( integer k = 1; k < n - 1; ++k )
    {
      J.insert( k, k - 1 ) = phi2_1( x( k - 1 ), x( k ) );
      J.insert( k, k )     = phi1_1( x( k ), x( k + 1 ) ) + phi2_2( x( k - 1 ), x( k ) );
      J.insert( k, k + 1 ) = phi1_2( x( k ), x( k + 1 ) );
    }
    J.insert( n - 1, n - 2 ) = phi2_1( x( n - 2 ), x( n - 1 ) );
    J.insert( n - 1, n - 1 ) = phi2_2( x( n - 2 ), x( n - 1 ) );
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.setZero();
  }
};
