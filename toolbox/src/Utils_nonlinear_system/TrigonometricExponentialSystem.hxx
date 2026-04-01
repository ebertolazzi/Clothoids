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

#define TRIGONOMETRIC_EXPONENTIAL_BIBTEX                                  \
  "Spedicato, E.\n"                                                       \
  "Computational experience with quasi-newton algoritms.\n"               \
  "for minimization problems of moderately large size.\n"                 \
  "Rep. CISE-N-175, Segrate (Milano), 1975.\n\n"                          \
  "@article{More:1981,\n"                                                 \
  "  author  = {Mor{\'e}, Jorge J. and Garbow, Burton S. and Hillstrom, " \
  "Kenneth E.},\n"                                                        \
  "  title   = {Testing Unconstrained Optimization Software},\n"          \
  "  journal = {ACM Trans. Math. Softw.},\n"                              \
  "  year    = {1981},\n"                                                 \
  "  volume  = {7},\n"                                                    \
  "  number  = {1},\n"                                                    \
  "  pages   = {17--41},\n"                                               \
  "  doi     = {10.1145/355934.355936},\n"                                \
  "}\n"

class TrigonometricExponentialSystem1 : public NonlinearSystem
{
public:
  TrigonometricExponentialSystem1( integer neq )
    : NonlinearSystem( "Trigonometric Exponential System prob 1", TRIGONOMETRIC_EXPONENTIAL_BIBTEX, neq )
  {
    check_even( n, 2 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    for ( integer i = 0; i < n; i += 2 )
      f( i ) = 3 * power3( x( i ) ) + 2 * x( i + 1 ) - 5 + sin( x( i ) - x( i + 1 ) ) * sin( x( i ) + x( i + 1 ) );
    for ( integer i = 1; i < n; i += 2 ) f( i ) = -x( i - 1 ) * exp( x( i - 1 ) - x( i ) ) + 4 * x( i ) - 3;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer i = 0; i < n; i += 2 )
    {
      J.insert( i, i )     = 9 * power2( x( i ) ) + sin( 2 * x( i ) );
      J.insert( i, i + 1 ) = 2 - sin( 2 * x( i + 1 ) );
    }
    for ( integer i = 1; i < n; i += 2 )
    {
      J.insert( i, i )     = x( i - 1 ) * exp( x( i - 1 ) - x( i ) ) + 4;
      J.insert( i, i - 1 ) = -( x( i - 1 ) + 1 ) * exp( x( i - 1 ) - x( i ) );
    }
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.setZero();
  }

  virtual void check_if_admissible( Vector const & x ) const override
  {
    for ( integer i = 0; i < n; ++i ) UTILS_ASSERT( std::abs( x( i ) ) < 100, "Bad range" );
  }

  virtual void bounding_box( Vector & L, Vector & U ) const override
  {
    U.fill( 100 );
    L.fill( -100 );
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class TrigonometricExponentialSystem2 : public NonlinearSystem
{
public:
  TrigonometricExponentialSystem2( integer neq )
    : NonlinearSystem( "Trigonometric Exponential System prob 2", TRIGONOMETRIC_EXPONENTIAL_BIBTEX, neq )
  {
    check_odd( n, 6 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = 3 * power3( x( 0 ) - x( 2 ) ) + 2 * x( 1 ) - 5 +
             sin( x( 0 ) - x( 1 ) - x( 2 ) ) * sin( x( 0 ) + x( 1 ) - x( 2 ) );

    for ( integer i = 1; i < n - 1; i += 2 )
      f( i ) = ( x( i + 1 ) - x( i - 1 ) ) * exp( x( i - 1 ) - x( i ) - x( i + 1 ) ) + 4 * x( i ) - 3;

    for ( integer i = 2; i < n - 1; i += 2 )
      f( i ) = 3 * power3( x( i ) - x( i + 2 ) ) + 6 * power3( x( i ) - x( i - 2 ) ) + 2 * x( i + 1 ) - 4 * x( i - 1 ) +
               5 - 2 * sin( x( i - 2 ) - x( i - 1 ) - x( i ) ) * sin( x( i - 2 ) + x( i - 1 ) - x( i ) ) +
               sin( x( i ) - x( i + 1 ) - x( i + 2 ) ) * sin( x( i ) + x( i + 1 ) - x( i + 2 ) );

    f( n - 1 ) = -6 * power3( x( n - 1 ) - x( n - 3 ) ) - 4 * x( n - 2 ) + 10 -
                 2 * sin( x( n - 3 ) - x( n - 2 ) - x( n - 1 ) ) * sin( x( n - 3 ) + x( n - 2 ) - x( n - 1 ) );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    J.insert( 0, 0 ) = 9 * power2( x( 0 ) - x( 2 ) ) + sin( 2 * ( x( 0 ) - x( 2 ) ) );
    J.insert( 0, 1 ) = 2 - sin( 2 * x( 1 ) );
    J.insert( 0, 2 ) = -9 * power2( x( 0 ) - x( 2 ) ) - sin( 2 * ( x( 0 ) - x( 2 ) ) );

    for ( integer i = 1; i < n - 1; i += 2 )
    {
      real_type ex         = exp( x( i - 1 ) - x( i ) - x( i + 1 ) );
      real_type tp         = x( i + 1 ) - x( i - 1 ) - 1;
      J.insert( i, i - 1 ) = tp * ex;
      J.insert( i, i )     = 4 + ( x( i - 1 ) - x( i + 1 ) ) * ex;
      J.insert( i, i + 1 ) = -tp * ex;
    }

    for ( integer i = 2; i < n - 1; i += 2 )
    {
      J.insert( i, i - 2 ) = 2 * sin( 2 * ( x( i ) - x( i - 2 ) ) ) - 18 * power2( x( i ) - x( i - 2 ) );
      J.insert( i, i - 1 ) = -4 + 2 * sin( 2 * x( i - 1 ) );
      J.insert( i, i - 0 ) = sin( 2 * ( x( i ) - x( i + 2 ) ) ) - 2 * sin( 2 * ( x( i ) - x( i - 2 ) ) ) +
                             27 * power2( x( i ) ) - 36 * x( i ) * x( i - 2 ) - 18 * x( i ) * x( i + 2 ) +
                             18 * power2( x( i - 2 ) ) + 9 * power2( x( i + 2 ) );
      J.insert( i, i + 1 ) = 2 - sin( 2 * x( i + 1 ) );
      J.insert( i, i + 2 ) = -9 * power2( x( i ) - x( i + 2 ) ) - sin( 2 * ( x( i ) - x( i + 2 ) ) );
    }

    J.insert( n - 1, n - 3 ) = -2 * sin( 2 * ( x( n - 3 ) - x( n - 1 ) ) ) + 18 * power2( x( n - 3 ) - x( n - 1 ) );
    J.insert( n - 1, n - 2 ) = 2 * sin( 2 * x( n - 2 ) ) - 4;
    J.insert( n - 1, n - 1 ) = 2 * sin( 2 * ( x( n - 3 ) - x( n - 1 ) ) ) - 18 * power2( x( n - 3 ) - x( n - 1 ) );
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 1 );
  }

  virtual void check_if_admissible( Vector const & x ) const override
  {
    for ( integer i = 0; i < n; ++i ) UTILS_ASSERT( std::abs( x( i ) ) < 100, "Bad range" );
  }

  virtual void bounding_box( Vector & L, Vector & U ) const override
  {
    U.fill( 100 );
    L.fill( -100 );
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class TrigExp : public NonlinearSystem
{
public:
  TrigExp( integer neq )
    : NonlinearSystem(
        "TrigExp",
        "@article{Ruggiero:1992,\n"
        "  author  = {Gomes-Ruggiero, M. and Martínez, J. and Moretti, "
        "A.},\n"
        "  title   = {Comparing Algorithms for Solving Sparse Nonlinear\n"
        "             Systems of Equations},\n"
        "  journal = {SIAM Journal on Scientific and Statistical "
        "Computing},\n"
        "  volume  = {13},\n"
        "  number  = {2},\n"
        "  pages   = {459-483},\n"
        "  year    = {1992},\n"
        "  doi     = {10.1137/0913025},\n"
        "}\n",
        neq )
  {
    check_min_equations( n, 3 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 )     = 3 * x( 0 ) * x( 0 ) + 2 * x( 1 ) - 5 + sin( x( 0 ) - x( 1 ) ) * sin( x( 0 ) + x( 1 ) );
    f( n - 1 ) = -x( n - 2 ) * exp( x( n - 1 ) - x( n - 2 ) ) + 4 * x( n - 1 ) - 3;
    for ( integer i = 1; i < n - 1; ++i )
      f( i ) = -x( i - 1 ) * exp( x( i - 1 ) - x( i ) ) + x( i ) * ( 4 + 3 * x( i ) * x( i ) ) + 2 * x( i + 1 ) +
               sin( x( i ) - x( i + 1 ) ) * sin( x( i ) + x( i + 1 ) );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    J.insert( 0, 0 ) = 6 * x( 0 ) + cos( x( 0 ) - x( 1 ) ) * sin( x( 0 ) + x( 1 ) ) +
                       sin( x( 0 ) - x( 1 ) ) * cos( x( 0 ) + x( 1 ) );
    J.insert( 0, 1 ) = 2 - cos( x( 0 ) - x( 1 ) ) * sin( x( 0 ) + x( 1 ) ) +
                       sin( x( 0 ) - x( 1 ) ) * cos( x( 0 ) + x( 1 ) );
    J.insert( n - 1, n - 1 ) = 4 - x( n - 2 ) * exp( x( n - 1 ) - x( n - 2 ) );
    J.insert( n - 1, n - 2 ) = ( x( n - 2 ) - 1 ) * exp( x( n - 1 ) - x( n - 2 ) );
    for ( integer i = 1; i < n - 1; ++i )
    {
      J.insert( i, i - 1 ) = -( 1 + x( i - 1 ) ) * exp( x( i - 1 ) - x( i ) );
      J.insert( i, i )     = x( i - 1 ) * exp( x( i - 1 ) - x( i ) ) + 9 * x( i ) * x( i ) + 4 +
                         cos( x( i ) - x( i + 1 ) ) * sin( x( i ) + x( i + 1 ) ) +
                         sin( x( i ) - x( i + 1 ) ) * cos( x( i ) + x( i + 1 ) );
      J.insert( i, i + 1 ) = 2 - cos( x( i ) - x( i + 1 ) ) * sin( x( i ) + x( i + 1 ) ) +
                             sin( x( i ) - x( i + 1 ) ) * cos( x( i ) + x( i + 1 ) );
    }
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 3 );
    auto & x0{ x_vec[0] };
    auto & x1{ x_vec[1] };
    auto & x2{ x_vec[2] };
    x0.resize( n );
    x1.resize( n );
    x2.resize( n );
    for ( integer k{ 0 }; k < n; ++k ) x0( k ) = 1 / real_type( n );
    x1.setZero();
    x2.fill( 0.3 );
  }

  virtual void check_if_admissible( Vector const & x ) const override
  {
    for ( integer i = 0; i < n; ++i ) UTILS_ASSERT( std::abs( x( i ) ) < 1000, "Bad range" );
  }

  virtual void bounding_box( Vector & L, Vector & U ) const override
  {
    U.fill( 1000 );
    L.fill( -1000 );
  }
};
