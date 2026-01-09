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

class GheriMancino : public NonlinearSystem
{
  real_type alpha, beta, gamma;

public:
  // sum log(xi-2)^2+log(xi-10)^2 - prod( xi) ^(1/5)
  GheriMancino( integer neq )
    : NonlinearSystem(
        "Gheri-Mancino function",
        "@Article{Gheri1971,\n"
        "  author  = {Gheri, G. and Mancino, O. G.},\n"
        "  title   = {A significant example to test methods for\n"
        "             solving systems of nonlinear equations},\n"
        "  journal = {CALCOLO},\n"
        "  year    = {1971},\n"
        "  volume  = {8},\n"
        "  number  = {1},\n"
        "  pages   = {107--113},\n"
        "  doi     = {10.1007/BF02575578}\n"
        "}\n",
        neq )
    , alpha( 7 )
    , beta( 17 )
    , gamma( 4 )
  {
  }

  real_type zfun( integer i, integer j, Vector const & x ) const
  {
    return sqrt( x( j ) * x( j ) + ( i + 1.0 ) / ( j + 1.0 ) );
  }

  real_type zfun_1( integer i, integer j, Vector const & x ) const
  {
    return x( j ) / sqrt( x( j ) * x( j ) + ( i + 1.0 ) / ( j + 1.0 ) );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    for ( integer i = 0; i < n; ++i )
    {
      f( i ) = beta * n * x( i ) + pow( i + 1 - 0.5 * n, gamma );
      for ( integer j = 0; j < n; ++j )
      {
        if ( i != j )
        {
          real_type zij = zfun( i, j, x );
          real_type lij = log( zij );
          real_type ss  = sin( lij );
          real_type cc  = cos( lij );
          f( i ) += zij * ( pow( ss, alpha ) + pow( cc, alpha ) );
        }
      }
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer i = 0; i < n; ++i )
    {
      for ( integer j = 0; j < n; ++j )
      {
        if ( i != j )
        {
          real_type zij    = zfun( i, j, x );
          real_type zij_1  = zfun_1( i, j, x );
          real_type lij    = log( zij );
          real_type ss     = sin( lij );
          real_type cc     = cos( lij );
          J.insert( i, j ) = zij_1 * ( pow( ss, alpha ) + pow( cc, alpha ) ) +
                             alpha * ( pow( ss, alpha - 1 ) * cc - pow( cc, alpha - 1 ) * ss ) * x( j ) / zij;
        }
        else
        {
          J.insert( i, j ) = beta * n;
        }
      }
    }
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    real_type c = beta * n - ( alpha + 1 ) * ( n - 1 );
    real_type k = beta * n + ( alpha + 1 ) * ( n - 1 );
    for ( integer i = 0; i < n; ++i )
    {
      x0( i ) = ( i + 0.5 * n ) * gamma;
      for ( integer j = 0; j < n; ++j )
      {
        if ( i != j )
        {
          real_type zij = sqrt( ( i + 1.0 ) / ( j + 1.0 ) );
          x0( i ) += zij * ( pow( sin( log( zij ) ), alpha ) + pow( cos( log( zij ) ), alpha ) );
        }
      }
      x0( i ) *= -( c + k ) / ( 2 * k );
    }
  }
};
