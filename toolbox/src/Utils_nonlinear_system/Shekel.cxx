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

#define SHEKEL_BIBTEX                                                  \
  "@book{brent2013,\n"                                                 \
  "  author    = {Brent, R.P.},\n"                                     \
  "  title     = {Algorithms for Minimization Without Derivatives},\n" \
  "  isbn      = {9780486143682},\n"                                   \
  "  series    = {Dover Books on Mathematics},\n"                      \
  "  year      = {2013},\n"                                            \
  "  publisher = {Dover Publications}\n"                               \
  "}\n"

class ShekelSQRN5 : public NonlinearSystem
{
  real_type c[5];
  real_type a[5][4];
  using Matrix = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;

public:
  // f = 0.0;
  // for j = 1 : m
  //   f = f - 1.0 / ( c(j) + sum ( ( x - a(1:n,j) ).^2 ) );
  // end

  ShekelSQRN5() : NonlinearSystem( "Shekel SQRN5 Function", SHEKEL_BIBTEX, 4 )
  {
    a[0][0] = 4.0;
    a[0][1] = 4.0;
    a[0][2] = 4.0;
    a[0][3] = 4.0;
    a[1][0] = 1.0;
    a[1][1] = 1.0;
    a[1][2] = 1.0;
    a[1][3] = 1.0;
    a[2][0] = 8.0;
    a[2][1] = 8.0;
    a[2][2] = 8.0;
    a[2][3] = 8.0;
    a[3][0] = 6.0;
    a[3][1] = 6.0;
    a[3][2] = 6.0;
    a[3][3] = 6.0;
    a[4][0] = 3.0;
    a[4][1] = 7.0;
    a[4][2] = 3.0;
    a[4][3] = 7.0;
    c[0]    = 0.1;
    c[1]    = 0.2;
    c[2]    = 0.2;
    c[3]    = 0.4;
    c[4]    = 0.6;
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = f( 1 ) = f( 2 ) = f( 3 ) = 0;
    for ( integer k = 0; k < n; ++k )
    {
      for ( integer j = 0; j < 5; ++j )
      {
        real_type d = c[j];
        for ( integer i = 0; i < n; ++i ) d += power2( x( i ) - a[i][j] );
        f( k ) += 2.0 * ( x( k ) - a[k][j] ) / ( d * d );
      }
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    Matrix J_full( n, n );
    J_full.setZero();
    for ( integer ii = 0; ii < n; ++ii )
    {
      for ( integer jj = 0; jj < n; ++jj )
      {
        for ( integer j = 0; j < 5; ++j )
        {
          real_type d = c[j];
          for ( integer i = 0; i < n; ++i ) d += power2( x( i ) - a[i][j] );
          J_full( ii, jj ) -= 8.0 * ( ( x[ii] - a[ii][j] ) * ( x[jj] - a[jj][j] ) ) / ( d * d * d );
          if ( ii == jj ) J_full( ii, jj ) += 2 / ( d * d );
        }
      }
    }
    J.resize( n, n );
    J = J_full.sparseView();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 4 );
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 1.0, 3.0, 5.0, 6.0;
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class ShekelSQRN7 : public NonlinearSystem
{
  real_type c[7];
  real_type a[7][4];
  using Matrix = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;

public:
  // f = 0.0;
  // for j = 1 : m
  //   f = f - 1.0 / ( c(j) + sum ( ( x - a(1:n,j) ).^2 ) );
  // end

  ShekelSQRN7() : NonlinearSystem( "Shekel SQRN7 Function", SHEKEL_BIBTEX, 4 )
  {
    a[0][0] = 4.0;
    a[0][1] = 4.0;
    a[0][2] = 4.0;
    a[0][3] = 4.0;
    a[1][0] = 1.0;
    a[1][1] = 1.0;
    a[1][2] = 1.0;
    a[1][3] = 1.0;
    a[2][0] = 8.0;
    a[2][1] = 8.0;
    a[2][2] = 8.0;
    a[2][3] = 8.0;
    a[3][0] = 6.0;
    a[3][1] = 6.0;
    a[3][2] = 6.0;
    a[3][3] = 6.0;
    a[4][0] = 3.0;
    a[4][1] = 7.0;
    a[4][2] = 3.0;
    a[4][3] = 7.0;
    a[5][0] = 2.0;
    a[5][1] = 9.0;
    a[5][2] = 2.0;
    a[5][3] = 9.0;
    a[6][0] = 5.0;
    a[6][1] = 5.0;
    a[6][2] = 3.0;
    a[6][3] = 3.0;
    c[0]    = 0.1;
    c[1]    = 0.2;
    c[2]    = 0.2;
    c[3]    = 0.4;
    c[4]    = 0.6;
    c[5]    = 0.6;
    c[6]    = 0.3;
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = f( 1 ) = f( 2 ) = f( 3 ) = 0;
    for ( integer k = 0; k < n; ++k )
    {
      for ( integer j = 0; j < 7; ++j )
      {
        real_type d = c[j];
        for ( integer i = 0; i < n; ++i ) d += power2( x( i ) - a[i][j] );
        f( k ) += 2.0 * ( x( k ) - a[k][j] ) / ( d * d );
      }
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    Matrix J_full( n, n );
    J_full.setZero();
    for ( integer ii = 0; ii < n; ++ii )
    {
      for ( integer jj = 0; jj < n; ++jj )
      {
        for ( integer j = 0; j < 7; ++j )
        {
          real_type d = c[j];
          for ( integer i = 0; i < n; ++i ) d += power2( x( i ) - a[i][j] );
          J_full( ii, jj ) -= 8.0 * ( ( x[ii] - a[ii][j] ) * ( x[jj] - a[jj][j] ) ) / ( d * d * d );
          if ( ii == jj ) J_full( ii, jj ) += 2 / ( d * d );
        }
      }
    }
    J.resize( n, n );
    J = J_full.sparseView();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 4 );
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 1.0, 3.0, 5.0, 6.0;
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class ShekelSQRN10 : public NonlinearSystem
{
  real_type c[10];
  real_type a[10][4];
  using Matrix = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;

public:
  // f = 0.0;
  // for j = 1 : m
  //   f = f - 1.0 / ( c(j) + sum ( ( x - a(1:n,j) ).^2 ) );
  // end

  ShekelSQRN10() : NonlinearSystem( "Shekel SQRN10 Function", SHEKEL_BIBTEX, 4 )
  {
    a[0][0] = 4.0;
    a[0][1] = 4.0;
    a[0][2] = 4.0;
    a[0][3] = 4.0;
    a[1][0] = 1.0;
    a[1][1] = 1.0;
    a[1][2] = 1.0;
    a[1][3] = 1.0;
    a[2][0] = 8.0;
    a[2][1] = 8.0;
    a[2][2] = 8.0;
    a[2][3] = 8.0;
    a[3][0] = 6.0;
    a[3][1] = 6.0;
    a[3][2] = 6.0;
    a[3][3] = 6.0;
    a[4][0] = 3.0;
    a[4][1] = 7.0;
    a[4][2] = 3.0;
    a[4][3] = 7.0;
    a[5][0] = 2.0;
    a[5][1] = 9.0;
    a[5][2] = 2.0;
    a[5][3] = 9.0;
    a[6][0] = 5.0;
    a[6][1] = 5.0;
    a[6][2] = 3.0;
    a[6][3] = 3.0;
    a[7][0] = 8.0;
    a[7][1] = 1.0;
    a[7][2] = 8.0;
    a[7][3] = 1.0;
    a[8][0] = 6.0;
    a[8][1] = 2.0;
    a[8][2] = 6.0;
    a[8][3] = 2.0;
    a[9][0] = 7.0;
    a[9][1] = 3.6;
    a[9][2] = 7.0;
    a[9][3] = 3.6;
    c[0]    = 0.1;
    c[1]    = 0.2;
    c[2]    = 0.2;
    c[3]    = 0.4;
    c[4]    = 0.6;
    c[5]    = 0.6;
    c[6]    = 0.3;
    c[7]    = 0.7;
    c[8]    = 0.5;
    c[9]    = 0.5;
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = f( 1 ) = f( 2 ) = f( 3 ) = 0;
    for ( integer k = 0; k < n; ++k )
    {
      for ( integer j = 0; j < 10; ++j )
      {
        real_type d = c[j];
        for ( integer i = 0; i < n; ++i ) d += power2( x( i ) - a[i][j] );
        f( k ) += 2.0 * ( x( k ) - a[k][j] ) / ( d * d );
      }
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    Matrix J_full( n, n );
    J_full.setZero();
    for ( integer ii = 0; ii < n; ++ii )
    {
      for ( integer jj = 0; jj < n; ++jj )
      {
        for ( integer j = 0; j < 10; ++j )
        {
          real_type d = c[j];
          for ( integer i = 0; i < n; ++i ) d += power2( x( i ) - a[i][j] );
          J_full( ii, jj ) -= 8.0 * ( ( x[ii] - a[ii][j] ) * ( x[jj] - a[jj][j] ) ) / ( d * d * d );
          if ( ii == jj ) J_full( ii, jj ) += 2 / ( d * d );
        }
      }
    }
    J.resize( n, n );
    J = J_full.sparseView();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 4 );
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 1.0, 3.0, 5.0, 6.0;
  }
};
