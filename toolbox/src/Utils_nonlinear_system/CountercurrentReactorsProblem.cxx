/*\
 |
 |  Author:
 |    Enrico Bertolazzi
 |    University of Trento
 |    Department of Industrial Engineering
 |    Via Sommarive 9, I-38123, Povo, Trento, Italy
 |    email: enrico.bertolazzi@unitn.it
\*/

#define COUNTERCURRENT_BIBTEX                                             \
  "@article{Bogle:1990,\n"                                                \
  "  author  = {Bogle, I. and Perkins, J.},\n"                            \
  "  title   = {A New Sparsity Preserving Quasi-Newton Update\n"          \
  "             for Solving Nonlinear Equations},\n"                      \
  "  journal = {SIAM Journal on Scientific and Statistical Computing},\n" \
  "  year    = {1990},\n"                                                 \
  "  volume  = {11},\n"                                                   \
  "  number  = {4},\n"                                                    \
  "  pages   = {621-630},\n"                                              \
  "  doi     = {10.1137/0911036},\n"                                      \
  "}\n"                                                                   \
  "@techreport{Bodon:1990,\n"                                             \
  "  author  = {Elena Bodon and Ladislav Luksan and Emilio Spedicato},\n" \
  "  title   = {Numerical performance of ABS codes for nonlinear least "  \
  "squares},\n"                                                           \
  "  year    = {2001},\n"                                                 \
  "  number  = {Tech. Rep. DMSIA 27/2001, Universita degli Studi di "     \
  "Bergamo}\n"                                                            \
  "}\n"

//  Problem N.2

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class CountercurrentReactorsProblem1 : public NonlinearSystem
{
  real_type const alpha;
  real_type const theta;

public:
  CountercurrentReactorsProblem1( integer neq )
    : NonlinearSystem( "Countercurrent Reactors Problem N.1", COUNTERCURRENT_BIBTEX, neq ), alpha( 0.5 ), theta( 4.0 )
  {
    check_min_equations( neq, 4 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const
  {
    for ( integer i = 0; i < n; i += 2 )
    {
      real_type xm2, xm1, xp2, xp3;
      if ( i == 0 )
      {
        xm2 = 1;
        xm1 = 0;
      }
      else
      {
        xm2 = x( i - 2 );
        xm1 = x( i - 1 );
      }
      if ( i >= n - 2 )
      {
        xp2 = 0;
        xp3 = 1;
      }
      else
      {
        xp2 = x( i + 2 );
        xp3 = x( i + 3 );
      }
      real_type xi  = x( i );
      real_type xp1 = x( i + 1 );
      f( i )        = alpha * xm2 + ( alpha - 1 ) * xp2 - xi * ( 1 + theta * xp1 );
      f( i + 1 )    = ( alpha - 1 ) * xm1 + ( alpha - 2 ) * xp3 - theta * xi * xp1;
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const
  {
    J.resize( n, n );
    J.setZero();
    for ( integer i = 0; i < n; i += 2 )
    {
      if ( i > 0 )
      {
        J.insert( i, i - 2 )     = alpha;
        J.insert( i + 1, i - 1 ) = alpha - 1;
      }
      if ( i < n - 3 )
      {
        J.insert( i, i + 2 )     = alpha - 1;
        J.insert( i + 1, i + 3 ) = alpha - 2;
      }
      real_type xi             = x( i );
      real_type xp1            = x( i + 1 );
      J.insert( i, i )         = -( 1 + theta * xp1 );
      J.insert( i, i + 1 )     = -theta * xi;
      J.insert( i + 1, i )     = -theta * xp1;
      J.insert( i + 1, i + 1 ) = -theta * xi;
    }
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const
  {
    x_vec.resize( 3 );
    x_vec[0].resize( n );
    x_vec[1].resize( n );
    x_vec[2].resize( n );

    auto & x0{ x_vec[0] };
    for ( integer i = 0; i < n; ++i )
    {
      switch ( i % 8 )
      {
        case 0: x0( i ) = 0.1; break;
        case 1:
        case 7: x0( i ) = 0.2; break;
        case 2:
        case 6: x0( i ) = 0.3; break;
        case 3:
        case 5: x0( i ) = 0.4; break;
        case 4: x0( i ) = 0.5; break;
      }
    }
    x_vec[1] = 10 * x0;
    x_vec[2] = 100 * x0;
  }
};

class CountercurrentReactorsProblem2 : public NonlinearSystem
{
  real_type const A0;
  real_type const A1;
  real_type const B0;
  real_type const theta;

public:
  CountercurrentReactorsProblem2( integer neq )
    : NonlinearSystem( "Countercurrent Reactors Problem N.2", COUNTERCURRENT_BIBTEX, neq )
    , A0( 1 )
    , A1( 0.414214 )
    , B0( 0 )
    , theta( 4 )
  {
    check_min_equations( neq, 6 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = A0 * x( 0 ) - ( 1 - x( 0 ) ) * x( 2 ) - A1 - theta * A1 * x( 1 );
    f( 1 ) = B0 * x( 0 ) - ( 1 - x( 0 ) ) * x( 3 ) - A1 - theta * A1 * x( 1 );
    f( 2 ) = A1 * x( 0 ) - ( 1 - x( 0 ) ) * x( 4 ) - x( 2 ) - theta * x( 2 ) * x( 3 );

    for ( integer i = 3; i < n; ++i )
    {
      real_type xp2 = 1;
      if ( i + 2 < n )
        xp2 = x( i + 2 );
      else if ( i + 2 == n )
        xp2 = 0;
      f( i ) = x( 0 ) * x( i - 2 ) - ( 1 - x( 0 ) ) * xp2 - x( i ) - theta * x( i - 1 ) * x( i );
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    J.insert( 0, 0 ) = A0 + x( 2 );
    J.insert( 0, 1 ) = -theta * A1;
    J.insert( 0, 2 ) = x( 0 ) - 1;

    J.insert( 1, 0 ) = B0 + x( 3 );
    J.insert( 1, 1 ) = -theta * A1;
    J.insert( 1, 3 ) = x( 0 ) - 1;

    J.insert( 2, 0 ) = A1 + x( 4 );
    J.insert( 2, 2 ) = -1 - theta * x( 3 );
    J.insert( 2, 3 ) = -theta * x( 2 );
    J.insert( 2, 4 ) = x( 0 ) - 1;

    for ( integer i = 3; i < n; ++i )
    {
      real_type xp2 = 1;
      if ( i + 2 < n )
        xp2 = x( i + 2 );
      else if ( i + 2 == n )
        xp2 = 0;
      J.insert( i, i - 2 ) = x( 0 );
      J.insert( i, i - 1 ) = -theta * x( i );
      J.insert( i, i )     = -1 - theta * x( i - 1 );
      if ( i + 2 < n ) J.insert( i, i + 2 ) = x( 0 ) - 1;  // Correzione: i+2 invece di i+1
      J.insert( i, 0 ) = x( i - 2 ) + xp2;
    }
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 3 );
    x_vec[0].resize( n );
    x_vec[1].resize( n );
    x_vec[2].resize( n );

    auto & x0{ x_vec[0] };
    for ( integer i = 0; i < n; ++i )
    {
      switch ( i % 8 )
      {
        case 0: x0( i ) = 0.1; break;
        case 1: x0( i ) = 0.2; break;
        case 2: x0( i ) = 0.3; break;
        case 3: x0( i ) = 0.4; break;
        case 4: x0( i ) = 0.5; break;
        case 5: x0( i ) = 0.4; break;
        case 6: x0( i ) = 0.3; break;
        case 7: x0( i ) = 0.2; break;
      }
    }
    x_vec[1] = 10 * x0;
    x_vec[2] = 100 * x0;
  }
};
