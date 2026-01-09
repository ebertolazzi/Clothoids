/*\
 |
 |  Author:
 |    Enrico Bertolazzi
 |    University of Trento
 |    Department of Industrial Engineering
 |    Via Sommarive 9, I-38123, Povo, Trento, Italy
 |    email: enrico.bertolazzi@unitn.it
\*/

#define PENALTY_FUNCTION_BIBTEX                                               \
  "@techreport{Raydan:2004,\n"                                                \
  "  author = {William La Cruz and Jose Mario Martinez and Marcos Raydan},\n" \
  "  title  = {Spectral residual method without gradient\n"                   \
  "             information for solving large-scale nonlinear\n"              \
  "             systems of equations: Theory and experiments},\n"             \
  "  number = {Technical Report RT-04-08},\n"                                 \
  "  year   = {2004}\n"                                                       \
  "}\n\n"                                                                     \
  "@article{LaCruz:2003,\n"                                                   \
  "  author    = {William {La Cruz}  and  Marcos Raydan},\n"                  \
  "  title     = {Nonmonotone Spectral Methods for Large-Scale Nonlinear "    \
  "Systems},\n"                                                               \
  "  journal   = {Optimization Methods and Software},\n"                      \
  "  year      = {2003},\n"                                                   \
  "  volume    = {18},\n"                                                     \
  "  number    = {5},\n"                                                      \
  "  pages     = {583--599},\n"                                               \
  "  publisher = {Taylor & Francis},\n"                                       \
  "  doi       = {10.1080/10556780310001610493},\n"                           \
  "}\n"

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class PenaltyIfunction : public NonlinearSystem
{
public:
  PenaltyIfunction( integer neq ) : NonlinearSystem( "Penalty I", PENALTY_FUNCTION_BIBTEX, neq )
  {
    check_min_equations( n, 2 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type sum = 0;
    for ( integer i = 0; i < n; ++i ) sum += x( i ) * x( i );
    for ( integer i = 0; i < n - 1; ++i ) f( i ) = sqrt( 1e-5 ) * ( x( i ) - 1 );
    f( n - 1 ) = ( sum / n - 1 ) / 4;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer i = 0; i < n - 1; ++i ) J.insert( i, i ) = sqrt( 1e-5 );
    for ( integer i = 0; i < n; ++i ) J.insert( n - 1, i ) = 0.5 * x( i ) / n;
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 1.0 / 3.0 );
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class PenaltyN1 : public NonlinearSystem
{
  real_type epsilon;

public:
  PenaltyN1( integer neq )
    : NonlinearSystem(
        "Penalty Function #1",
        "@book{brent2013,\n"
        "  author    = {Brent, R.P.},\n"
        "  title     = {Algorithms for Minimization Without Derivatives},\n"
        "  isbn      = {9780486143682},\n"
        "  series    = {Dover Books on Mathematics},\n"
        "  year      = {2013},\n"
        "  publisher = {Dover Publications}\n"
        "}\n",
        neq )
    , epsilon( 0.00001 )
  {
    check_min_equations( n, 2 );
  }

  real_type sum( Vector const & x ) const
  {
    real_type t1 = 0;
    for ( integer i = 0; i < n; ++i ) t1 += x( i ) * x( i );
    return 4 * t1 - 1;
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type ap = 2 * epsilon;
    real_type t1 = sum( x );
    for ( integer i = 0; i < n; ++i ) f( i ) = ( ap + t1 ) * x( i ) - ap;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type ap = 2 * epsilon;
    real_type t1 = sum( x );
    for ( integer i = 0; i < n; ++i )
    {
      for ( integer j = 0; j < n; ++j )
      {
        real_type tmp = 8 * x( j ) * x( i );
        if ( i == j ) tmp += ap + t1;
        J.insert( i, j ) = tmp;
      }
    }
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    for ( integer i{ 0 }; i < n; ++i ) x0( i ) = i + 1;
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class PenaltyN2 : public NonlinearSystem
{
  using Matrix = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;

  real_type epsilon;

public:
  PenaltyN2( integer neq )
    : NonlinearSystem(
        "Penalty Function #2",
        "@book{brent2013,\n"
        "  author    = {Brent, R.P.},\n"
        "  title     = {Algorithms for Minimization Without Derivatives},\n"
        "  isbn      = {9780486143682},\n"
        "  series    = {Dover Books on Mathematics},\n"
        "  year      = {2013},\n"
        "  publisher = {Dover Publications}\n"
        "}\n",
        neq )
    , epsilon( 0.00001 )
  {
    check_min_equations( n, 2 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type ap = epsilon;

    real_type t1 = -1.0;
    for ( integer j = 0; j < n; ++j ) t1 += ( n - j ) * ( x( j ) * x( j ) );

    real_type d2 = 1.0;
    real_type th = 4.0 * t1;
    real_type s2 = 0.0;
    for ( integer j = 0; j < n; ++j )
    {
      f( j )       = ( n - j ) * x( j ) * th;
      real_type s1 = exp( x( j ) / 10.0 );
      if ( j > 0 )
      {
        real_type s3 = s1 + s2 - d2 * ( exp( 0.1 ) + 1.0 );
        f( j ) += ap * s1 * ( s3 + s1 - 1.0 / exp( 0.1 ) ) / 5.0;
        f( j - 1 ) += ap * s2 * s3 / 5.0;
      }
      s2 = s1;
      d2 = d2 * exp( 0.1 );
    }
    f( 0 ) += 2.0 * ( x( 0 ) - 0.2 );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    Matrix J_full( n, n );
    J_full.setZero();

    real_type ap = 2 * epsilon;

    real_type t1 = -1.0;
    for ( integer j = 0; j < n; ++j ) t1 += ( n - j ) * ( x( j ) * x( j ) );

    real_type d1 = exp( 0.1 );
    real_type d2 = 1.0;
    real_type s2 = 0.0;
    real_type th = 4.0 * t1;

    for ( integer j = 0; j < n; ++j )
    {
      J_full( j, j ) = 8.0 * power2( ( n - j ) * x( j ) ) + ( n - j ) * th;

      real_type s1 = exp( x( j ) / 10.0 );

      if ( j > 0 )
      {
        real_type s3 = s1 + s2 - d2 * ( d1 + 1.0 );
        J_full( j, j ) += ap * s1 * ( s3 + s1 - 1.0 / d1 + 2.0 * s1 ) / 50.0;
        J_full( j - 1, j - 1 ) += ap * s2 * ( s2 + s3 ) / 50.0;
        for ( integer k = 0; k < j; ++k ) { J_full( j, k ) = 8.0 * ( n - j ) * ( n - k ) * x( j ) * x( k ); }
        J_full( j, j - 1 ) += ap * s1 * s2 / 50.0;
      }
      s2 = s1;
      d2 = d1 * d2;
    }

    J_full( 0, 0 ) += 2;
    for ( integer i = 0; i < n; ++i )
      for ( integer j = i + 1; j < n; ++j ) J_full( i, j ) = J_full( j, i );
    J.resize( n, n );
    J = J_full.sparseView();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 0.5 );
  }
};
