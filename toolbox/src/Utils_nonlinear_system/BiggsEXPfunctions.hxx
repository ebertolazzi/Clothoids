/*\
 |
 |  Author:
 |    Enrico Bertolazzi
 |    University of Trento
 |    Department of Industrial Engineering
 |    Via Sommarive 9, I-38123, Povo, Trento, Italy
 |    email: enrico.bertolazzi@unitn.it
\*/

/*
  (1/2) sum (eq(k,x)-yk)^2

  -->

  sum (eq(k,x)-yk) * D_x(j) eq(k,x)

*/

/*
  (1/2) sum (e(k)-yk)^2 + sum l(k)*(e(k)-eq(k,x))

  e(k)-zk+l(k) = 0
  e(k)-eq(k,x) = 0
  -sum_k D_x(j) eq(k,x) * l(k) = 0

  // si può eliminare l(k)

  eq(k,x)-e(k) = 0
  sum_k D_x(j) eq(k,x) * (e(k)-yk) = 0    j = 1,2,.,,nx
*/

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

#define BIGGS_BIBTEX                                                      \
  "@Article{Biggs:1971,\n"                                                \
  "  author = {M. C. Biggs},\n"                                           \
  "  title  = { Minimization algorithms making use of non-quadratic\n"    \
  "             properties of the objective function},\n"                 \
  "  volume = {8},\n"                                                     \
  "  pages  = {315--327},\n"                                              \
  "  year   = {1971},\n"                                                  \
  "  journal = {Journal of the Institute of Mathematics and its "         \
  "Applications}\n"                                                       \
  "}\n\n"                                                                 \
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

class BiggsEXP2function : public NonlinearSystem
{
  Vector        z, y;
  integer const NPT;

  using Matrix = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;

public:
  BiggsEXP2function() : NonlinearSystem( "Biggs EXP2 function", BIGGS_BIBTEX, 2 ), NPT( 10 )
  {
    z.resize( NPT );
    y.resize( NPT );
    for ( integer i = 0; i < NPT; ++i ) z( i ) = ( i + 1 ) * 0.1;
    y = exp( -z.array() ) - 5 * exp( -10 * z.array() );
  }

  void map( Vector const & x, Vector & eq ) const
  {
    eq = exp( -x( 0 ) * z.array() ) - 5 * exp( -x( 1 ) * z.array() ) - y.array();
  }

  void Grad_map( Vector const & x, integer k, Vector & G ) const
  {
    G( 0 ) = -z( k ) * exp( -x( 0 ) * z( k ) );
    G( 1 ) = 5 * z( k ) * exp( -x( 1 ) * z( k ) );
  }

  void Hess_map( Vector const & x, integer k, Matrix & H ) const
  {
    real_type zk  = z( k );
    real_type zk2 = zk * zk;
    real_type ex0 = exp( -x( 0 ) * zk );
    real_type ex1 = exp( -x( 1 ) * zk );
    H( 0, 0 )     = zk2 * ex0;
    H( 0, 1 )     = 0;
    H( 1, 1 )     = -5 * zk2 * ex1;
    H( 1, 0 )     = 0;
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    Vector eq( NPT ), G( n );
    map( x, eq );
    f.setZero();
    for ( integer k = 0; k < NPT; ++k )
    {
      Grad_map( x, k, G );
      f += eq( k ) * G;
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    Vector eq( NPT ), G( n );
    Matrix H( n, n ), J_full( n, n );
    map( x, eq );
    J_full.setZero();
    for ( integer k = 0; k < NPT; ++k )
    {
      Grad_map( x, k, G );
      Hess_map( x, k, H );
      for ( integer i = 0; i < n; ++i )
      {
        for ( integer j = 0; j < n; ++j ) { J_full( i, j ) += eq( k ) * H( i, j ) + G( i ) * G( j ); }
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
    x0 << 16.7046761257, 16.7046761257;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 1, 2;
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class BiggsEXP3function : public NonlinearSystem
{
  Vector        z, y;
  integer const NPT;

  using Matrix = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;

public:
  BiggsEXP3function() : NonlinearSystem( "Biggs EXP3 function", BIGGS_BIBTEX, 3 ), NPT( 10 )
  {
    z.resize( NPT );
    y.resize( NPT );
    for ( integer i = 0; i < NPT; ++i ) z( i ) = ( i + 1 ) * 0.1;
    y = exp( -z.array() ) - 5 * exp( -10 * z.array() );
  }

  void map( Vector const & x, Vector & eq ) const
  {
    eq = exp( -x( 0 ) * z.array() ) - x( 2 ) * exp( -x( 1 ) * z.array() ) - y.array();
  }

  void Grad_map( Vector const & x, integer k, Vector & G ) const
  {
    real_type zk  = z( k );
    real_type ex0 = exp( -x( 0 ) * zk );
    real_type ex1 = exp( -x( 1 ) * zk );
    G( 0 )        = -zk * ex0;
    G( 1 )        = x( 2 ) * zk * ex1;
    G( 2 )        = -ex1;
  }

  void Hess_map( Vector const & x, integer k, Matrix & H ) const
  {
    real_type zk  = z( k );
    real_type zk2 = zk * zk;
    real_type ex0 = exp( -x( 0 ) * zk );
    real_type ex1 = exp( -x( 1 ) * zk );
    H( 0, 0 )     = zk2 * ex0;
    H( 0, 1 )     = 0;
    H( 0, 2 )     = 0;
    H( 1, 0 )     = 0;
    H( 1, 1 )     = -x( 2 ) * zk2 * ex1;
    H( 1, 2 )     = zk * ex1;
    H( 2, 0 )     = 0;
    H( 2, 1 )     = zk * ex1;
    H( 2, 2 )     = 0;
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    Vector eq( NPT ), G( n );
    map( x, eq );
    f.setZero();
    for ( integer k = 0; k < NPT; ++k )
    {
      Grad_map( x, k, G );
      f += eq( k ) * G;
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    Vector eq( NPT ), G( n );
    Matrix H( n, n ), J_full( n, n );
    map( x, eq );
    J_full.setZero();
    for ( integer k = 0; k < NPT; ++k )
    {
      Grad_map( x, k, G );
      Hess_map( x, k, H );
      for ( integer i = 0; i < n; ++i )
      {
        for ( integer j = 0; j < n; ++j ) { J_full( i, j ) += eq( k ) * H( i, j ) + G( i ) * G( j ); }
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
    x0 << 1, 1, 5;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 1, 2, 1;
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class BiggsEXP4function : public NonlinearSystem
{
  Vector        z, y;
  integer const NPT;

  using Matrix = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;

public:
  BiggsEXP4function() : NonlinearSystem( "Biggs EXP4 function", BIGGS_BIBTEX, 4 ), NPT( 10 )
  {
    z.resize( NPT );
    y.resize( NPT );
    for ( integer i = 0; i < NPT; ++i ) z( i ) = ( i + 1 ) * 0.1;
    y = exp( -z.array() ) - 5 * exp( -10 * z.array() );
  }

  void map( Vector const & x, Vector & eq ) const
  {
    eq = x( 2 ) * exp( -x( 0 ) * z.array() ) - x( 3 ) * exp( -x( 1 ) * z.array() ) - y.array();
  }

  void Grad_map( Vector const & x, integer k, Vector & G ) const
  {
    real_type zk  = z( k );
    real_type ex0 = exp( -x( 0 ) * zk );
    real_type ex1 = exp( -x( 1 ) * zk );
    G( 0 )        = -x( 2 ) * zk * ex0;
    G( 1 )        = x( 3 ) * zk * ex1;
    G( 2 )        = ex0;
    G( 3 )        = -ex1;
  }

  void Hess_map( Vector const & x, integer k, Matrix & H ) const
  {
    real_type zk  = z( k );
    real_type zk2 = zk * zk;
    real_type ex0 = exp( -x( 0 ) * zk );
    real_type ex1 = exp( -x( 1 ) * zk );

    H( 0, 0 ) = x( 2 ) * zk2 * ex0;
    H( 0, 1 ) = 0;
    H( 0, 2 ) = -zk * ex0;
    H( 0, 3 ) = 0;

    H( 1, 0 ) = 0;
    H( 1, 1 ) = -x( 3 ) * zk2 * ex1;
    H( 1, 2 ) = 0;
    H( 1, 3 ) = zk * ex1;

    H( 2, 0 ) = -zk * ex0;
    H( 2, 1 ) = 0;
    H( 2, 2 ) = 0;
    H( 2, 3 ) = 0;

    H( 3, 0 ) = 0;
    H( 3, 1 ) = zk * ex1;
    H( 3, 2 ) = 0;
    H( 3, 3 ) = 0;
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    Vector eq( NPT ), G( n );
    map( x, eq );
    f.setZero();
    for ( integer k = 0; k < NPT; ++k )
    {
      Grad_map( x, k, G );
      f += eq( k ) * G;
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    Vector eq( NPT ), G( n );
    Matrix H( n, n ), J_full( n, n );
    map( x, eq );
    J_full.setZero();
    for ( integer k = 0; k < NPT; ++k )
    {
      Grad_map( x, k, G );
      Hess_map( x, k, H );
      for ( integer i = 0; i < n; ++i )
      {
        for ( integer j = 0; j < n; ++j ) { J_full( i, j ) += eq( k ) * H( i, j ) + G( i ) * G( j ); }
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
    x0 << -1.370515321, -1.370515321, 0.146554089168054, 0.003045452799256;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 1, 2, 1, 1;
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class BiggsEXP5function : public NonlinearSystem
{
  Vector        z, y;
  integer const NPT;

  using Matrix = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;

public:
  BiggsEXP5function() : NonlinearSystem( "Biggs EXP5 function", BIGGS_BIBTEX, 5 ), NPT( 11 )
  {
    z.resize( NPT );
    y.resize( NPT );
    for ( integer k = 0; k < NPT; ++k ) z( k ) = ( k + 1 ) * 0.1;
    y = exp( -z.array() ) - 5 * exp( -10 * z.array() ) + 3 * exp( -4 * z.array() );
  }

  void map( Vector const & x, Vector & eq ) const
  {
    eq = x( 2 ) * exp( -x( 0 ) * z.array() ) - x( 3 ) * exp( -x( 1 ) * z.array() ) + 3 * exp( -x( 4 ) * z.array() ) -
         y.array();
  }

  void Grad_map( Vector const & x, integer k, Vector & G ) const
  {
    real_type zk  = z( k );
    real_type ex0 = exp( -x( 0 ) * zk );
    real_type ex1 = exp( -x( 1 ) * zk );
    real_type ex4 = exp( -x( 4 ) * zk );
    G( 0 )        = -x( 2 ) * zk * ex0;
    G( 1 )        = x( 3 ) * zk * ex1;
    G( 2 )        = ex0;
    G( 3 )        = -ex1;
    G( 4 )        = -3 * zk * ex4;
  }

  void Hess_map( Vector const & x, integer k, Matrix & H ) const
  {
    real_type zk  = z( k );
    real_type zk2 = zk * zk;
    real_type ex0 = exp( -x( 0 ) * zk );
    real_type ex1 = exp( -x( 1 ) * zk );
    real_type ex4 = exp( -x( 4 ) * zk );

    H( 0, 0 ) = x( 2 ) * zk2 * ex0;
    H( 0, 1 ) = 0;
    H( 0, 2 ) = -zk * ex0;
    H( 0, 3 ) = 0;
    H( 0, 4 ) = 0;

    H( 1, 0 ) = 0;
    H( 1, 1 ) = -x( 3 ) * zk2 * ex1;
    H( 1, 2 ) = 0;
    H( 1, 3 ) = zk * ex1;
    H( 1, 4 ) = 0;

    H( 2, 0 ) = -zk * ex0;
    H( 2, 1 ) = 0;
    H( 2, 2 ) = 0;
    H( 2, 3 ) = 0;
    H( 2, 4 ) = 0;

    H( 3, 0 ) = 0;
    H( 3, 1 ) = zk * ex1;
    H( 3, 2 ) = 0;
    H( 3, 3 ) = 0;
    H( 3, 4 ) = 0;

    H( 4, 0 ) = 0;
    H( 4, 1 ) = 0;
    H( 4, 2 ) = 0;
    H( 4, 3 ) = 0;
    H( 4, 4 ) = 3 * zk2 * ex4;
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    Vector eq( NPT ), G( n );
    map( x, eq );
    f.setZero();
    for ( integer k = 0; k < NPT; ++k )
    {
      Grad_map( x, k, G );
      f += eq( k ) * G;
    }
  }

  // sum (eq(k,x)-yk) * D_x(j) eq(k,x)
  // sum (eq(k,x)-yk)^2 - z

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    Vector eq( NPT ), G( n );
    Matrix H( n, n ), J_full( n, n );
    map( x, eq );
    J_full.setZero();
    for ( integer k = 0; k < NPT; ++k )
    {
      Grad_map( x, k, G );
      Hess_map( x, k, H );
      for ( integer i = 0; i < n; ++i )
      {
        for ( integer j = 0; j < n; ++j ) { J_full( i, j ) += eq( k ) * H( i, j ) + G( i ) * G( j ); }
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
    x0 << 2.645310671110160, 2.693935021162610, 0.048826924681352, 0.049142366883430, 0.025390267672614;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 1, 2, 1, 1, 1;
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class BiggsEXP6function : public NonlinearSystem
{
  real_type const xmin;
  real_type const xmax;
  Vector          y, t;

public:
  BiggsEXP6function() : NonlinearSystem( "Biggs EXP6 function", BIGGS_BIBTEX, 6 ), xmin( -10 ), xmax( 20 )
  {
    y.resize( 6 );
    t.resize( 6 );
    for ( integer i = 0; i < 6; ++i ) t( i ) = ( i + 1 ) * 0.1;
    y = exp( -t.array() ) - 5 * exp( -10 * t.array() ) + 3 * exp( -4 * t.array() );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    for ( integer k = 0; k < 6; ++k )
    {
      real_type e0 = exp( -t( k ) * x( 0 ) );
      real_type e1 = exp( -t( k ) * x( 1 ) );
      real_type e4 = exp( -t( k ) * x( 4 ) );
      f( k )       = x( 2 ) * e0 - x( 3 ) * e1 + x( 5 ) * e4 - y( k );
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer k = 0; k < 6; ++k )
    {
      real_type e0     = exp( -t( k ) * x( 0 ) );
      real_type e1     = exp( -t( k ) * x( 1 ) );
      real_type e4     = exp( -t( k ) * x( 4 ) );
      J.insert( k, 0 ) = -x( 2 ) * t( k ) * e0;
      J.insert( k, 1 ) = x( 3 ) * t( k ) * e1;
      J.insert( k, 2 ) = e0;
      J.insert( k, 3 ) = -e1;
      J.insert( k, 4 ) = -x( 5 ) * t( k ) * e4;
      J.insert( k, 5 ) = e4;
    }
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 2 );
    auto & x0{ x_vec[0] };
    auto & x1{ x_vec[1] };
    x0.resize( n );
    x1.resize( n );

    x0 << 10, 4, -5, -3, 1, 1;
    x1 << 1, 10, 1, 5, 4, 3;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 1, 2, 1, 1, 1, 1;
  }

  virtual void check_if_admissible( Vector const & x ) const override
  {
    UTILS_ASSERT(
      x( 0 ) > xmin && x( 0 ) < xmax && x( 1 ) > xmin && x( 1 ) < xmax && x( 2 ) > xmin && x( 2 ) < xmax &&
        x( 3 ) > xmin && x( 3 ) < xmax && x( 4 ) > xmin && x( 4 ) < xmax && x( 5 ) > xmin && x( 5 ) < xmax,
      "Bad Range" );
  }

  virtual void bounding_box( Vector & L, Vector & U ) const override
  {
    L.fill( xmin );
    U.fill( xmax );
  }
};
