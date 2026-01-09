/*\
 |
 |  Author:
 |    Enrico Bertolazzi
 |    University of Trento
 |    Department of Industrial Engineering
 |    Via Sommarive 9, I-38123, Povo, Trento, Italy
 |    email: enrico.bertolazzi@unitn.it
\*/

#define DEVILLIERS_BIBTEX                                                 \
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

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class DeVilliersGlasser01 : public NonlinearSystem
{
  real_type t_vec[24], y_vec[24];

  using Matrix = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;

public:
  DeVilliersGlasser01() : NonlinearSystem( "De Villiers Glasser Problem #1", DEVILLIERS_BIBTEX, 4 )
  {
    for ( integer i = 0; i < 24; ++i )
    {
      t_vec[i] = i / 10.0;
      y_vec[i] = 60.137 * pow( 1.371, t_vec[i] ) * sin( 3.112 * t_vec[i] + 1.761 );
    }
  }

  real_type f_val( Vector const & x, integer k ) const
  {
    real_type t  = t_vec[k];
    real_type bf = x( 0 ) * pow( x( 1 ), t );
    return sin( x( 2 ) * t + x( 3 ) ) * bf - y_vec[k];
  }

  void f_grad( Vector const & x, integer k, Vector & g ) const
  {
    real_type t   = t_vec[k];
    real_type SIN = sin( x( 2 ) * t + x( 3 ) );
    real_type COS = cos( x( 2 ) * t + x( 3 ) );
    real_type PW  = pow( x( 1 ), t );
    real_type PW1 = pow( x( 1 ), t - 1 ) * t;
    g( 0 )        = SIN * PW;
    g( 1 )        = SIN * x( 0 ) * PW1;
    g( 2 )        = COS * t * x( 0 ) * PW;
    g( 3 )        = COS * x( 0 ) * PW;
  }

  void f_hess( Vector const & x, integer k, Matrix & h ) const
  {
    real_type t   = t_vec[k];
    real_type SIN = sin( x( 2 ) * t + x( 3 ) );
    real_type COS = cos( x( 2 ) * t + x( 3 ) );
    real_type PW  = pow( x( 1 ), t );
    real_type PW1 = pow( x( 1 ), t - 1 ) * t;
    real_type PW2 = pow( x( 1 ), t - 2 ) * t * ( t - 1 );
    h( 0, 0 )     = 0;
    h( 0, 1 )     = SIN * PW1;
    h( 0, 2 )     = COS * PW * t;
    h( 0, 3 )     = COS * PW;

    h( 1, 0 ) = h( 0, 1 );
    h( 1, 1 ) = SIN * x( 0 ) * PW2;
    h( 1, 2 ) = COS * x( 0 ) * PW1 * t;
    h( 1, 3 ) = COS * x( 0 ) * PW1;

    h( 2, 0 ) = h( 0, 2 );
    h( 2, 1 ) = h( 1, 2 );
    h( 2, 2 ) = -SIN * t * t * x( 0 ) * PW;
    h( 2, 3 ) = -SIN * t * x( 0 ) * PW;

    h( 3, 0 ) = h( 0, 3 );
    h( 3, 1 ) = h( 1, 3 );
    h( 3, 2 ) = h( 2, 3 );
    h( 3, 3 ) = -SIN * x( 0 ) * PW;
  }

  bool check_x( Vector const & x ) const { return x( 1 ) > 0; }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    if ( check_x( x ) )
    {
      Vector g( 4 );
      f.setZero();
      for ( integer k = 0; k < 24; ++k )
      {
        real_type ff = f_val( x, k );
        f_grad( x, k, g );
        f += ff * g;
      }
    }
    else
    {
      f.fill( nan( "DeVilliersGlasser01" ) );
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    Vector g( 4 );
    Matrix h( 4, 4 ), J_full( 4, 4 );
    J_full.setZero();
    for ( integer k = 0; k < 24; ++k )
    {
      real_type ff = f_val( x, k );
      f_grad( x, k, g );
      f_hess( x, k, h );
      for ( integer ii = 0; ii < 4; ++ii )
        for ( integer jj = 0; jj < 4; ++jj ) J_full( ii, jj ) += ff * h( ii, jj ) + g( ii ) * g( jj );
    }
    J.resize( n, n );
    J = J_full.sparseView();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 1, 8, 4, 4.412;
  }

  virtual void check_if_admissible( Vector const & x ) const override
  {
    UTILS_ASSERT( x( 1 ) > 0, "DeVilliersGlasser01, x(1) = {} <= 0", x( 1 ) );
    UTILS_ASSERT( -500 <= x( 0 ) && x( 0 ) <= 500, "DeVilliersGlasser#1, x(0) = {} must be in [-500,500]", x( 0 ) );
    UTILS_ASSERT( 0 <= x( 1 ) && x( 1 ) <= 500, "DeVilliersGlasser#1, x(1) = {} must be in [0,500]", x( 1 ) );
    UTILS_ASSERT( -500 <= x( 2 ) && x( 2 ) <= 500, "DeVilliersGlasser#1, x(2) = {} must be in [-500,500]", x( 2 ) );
    UTILS_ASSERT( -500 <= x( 3 ) && x( 3 ) <= 500, "DeVilliersGlasser#1, x(3) = {} must be in [-500,500]", x( 3 ) );
  }

  virtual void bounding_box( Vector & L, Vector & U ) const override
  {
    L[0] = -500;
    U[0] = 500;
    L[1] = 0;
    U[1] = 500;
    L[2] = -500;
    U[2] = 500;
    L[3] = -500;
    U[3] = 500;
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class DeVilliersGlasser02 : public NonlinearSystem
{
  real_type t_vec[24], y_vec[24];

  using Matrix = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;

public:
  DeVilliersGlasser02() : NonlinearSystem( "De Villiers Glasser Problem #2", DEVILLIERS_BIBTEX, 5 )
  {
    for ( integer i = 0; i < 24; ++i )
    {
      t_vec[i] = i / 10.0;
      y_vec[i] = 53.81 * pow( 1.27, t_vec[i] ) * tanh( 3.012 * t_vec[i] + sin( 2.13 * t_vec[i] ) ) *
                 cos( exp( 0.507 ) * t_vec[i] );
    }
  }

  real_type f_val( Vector const & x, integer k ) const
  {
    real_type t = t_vec[k];
    real_type bf;
    if ( x( 1 ) < 0 )
      bf = pow( -x( 1 ), t );
    else
      bf = pow( x( 1 ), t );
    return cos( exp( x( 4 ) ) * t ) * tanh( x( 2 ) * t + sin( x( 3 ) * t ) ) * bf * x( 0 ) - y_vec[k];
  }

  void f_grad( Vector const & x, integer k, Vector & g ) const
  {
    real_type PW   = pow( std::abs( x( 1 ) ), t_vec[k] );
    real_type PW1  = pow( std::abs( x( 1 ) ), t_vec[k] - 1 );
    real_type TANH = tanh( x( 2 ) * t_vec[k] + sin( x( 3 ) * t_vec[k] ) );
    real_type COSH = cosh( x( 2 ) * t_vec[k] + sin( x( 3 ) * t_vec[k] ) );
    real_type COS  = cos( t_vec[k] * exp( x( 4 ) ) );
    real_type SIN  = sin( t_vec[k] * exp( x( 4 ) ) );
    real_type bf   = x( 0 ) * t_vec[k];
    if ( x( 1 ) < 0 ) PW1 = -PW1;
    g( 0 ) = PW * TANH * COS;
    g( 1 ) = bf * PW1 * TANH * COS;
    g( 2 ) = bf * PW * COS / ( COSH * COSH );
    g( 3 ) = bf * PW * cos( x( 3 ) * t_vec[k] ) * COS / ( COSH * COSH );
    g[4]   = -bf * PW * TANH * exp( x( 4 ) ) * SIN;
  }

  void f_hess( Vector const & x, integer k, Matrix & h ) const
  {
    real_type PW   = pow( std::abs( x( 1 ) ), t_vec[k] );
    real_type PW1  = pow( std::abs( x( 1 ) ), t_vec[k] - 1 );
    real_type PW2  = pow( std::abs( x( 1 ) ), t_vec[k] - 2 );
    real_type TANH = tanh( x( 2 ) * t_vec[k] + sin( x( 3 ) * t_vec[k] ) );
    real_type COSH = cosh( x( 2 ) * t_vec[k] + sin( x( 3 ) * t_vec[k] ) );
    real_type SINH = sinh( x( 2 ) * t_vec[k] + sin( x( 3 ) * t_vec[k] ) );
    real_type COS  = cos( t_vec[k] * exp( x( 4 ) ) );
    real_type SIN  = sin( t_vec[k] * exp( x( 4 ) ) );
    if ( x( 1 ) < 0 ) PW1 = -PW1;
    if ( x( 1 ) < 0 ) PW2 = -PW2;
    real_type t1  = t_vec[k];
    real_type t4  = PW1 * t1 * TANH * COS;
    real_type t5  = t1 * PW;
    real_type t6  = COSH * COSH;
    real_type t7  = 1 / t6;
    real_type t9  = t5 * COS * t7;
    real_type t10 = x( 3 ) * t1;
    real_type t11 = cos( t10 );
    real_type t13 = t11 * COS * t7;
    real_type t14 = t5 * t13;
    real_type t16 = exp( x( 4 ) );
    real_type t17 = t1 * t16;
    real_type t19 = PW * TANH * t17 * SIN;
    real_type t23 = SINH / COSH;
    real_type t28 = x( 0 ) * PW1;
    real_type t29 = t1 * t1;
    real_type t32 = t7 * COS * t29 * t28;
    real_type t33 = t29 * t28;
    real_type t34 = t13 * t33;
    real_type t37 = SIN * t16 * TANH * t33;
    real_type t39 = t29 * PW * x( 0 );
    real_type t42 = 1 / t6 / COSH;
    real_type t50 = 2.0 * t42 * COS * t11 * SINH * t39;
    real_type t53 = t7 * SIN * t16 * t39;
    real_type t54 = t11 * t11;
    real_type t57 = sin( t10 );
    real_type t68 = t7 * SIN * t11 * t16 * t39;
    h( 0, 0 )     = 0.0;
    h( 0, 1 )     = t4;
    h( 0, 2 )     = t9;
    h( 0, 3 )     = t14;
    h( 0, 4 )     = -t19;
    h( 1, 0 )     = t4;
    h( 1, 1 )     = ( t1 - 1.0 ) * COS * t23 * x( 0 ) * PW2 * t1;
    h( 1, 2 )     = t32;
    h( 1, 3 )     = t34;
    h( 1, 4 )     = -t37;
    h( 2, 0 )     = t9;
    h( 2, 1 )     = t32;
    h( 2, 2 )     = -2.0 * t42 * COS * SINH * t39;
    h( 2, 3 )     = -t50;
    h( 2, 4 )     = -t53;
    h( 3, 0 )     = t14;
    h( 3, 1 )     = t34;
    h( 3, 2 )     = -t50;
    h( 3, 3 )     = -t42 * t29 * PW * x( 0 ) * ( COSH * t57 + 2.0 * SINH * t54 ) * COS;
    h( 3, 4 )     = -t68;
    h( 4, 0 )     = -t19;
    h( 4, 1 )     = -t37;
    h( 4, 2 )     = -t53;
    h( 4, 3 )     = -t68;
    h( 4, 4 )     = -( COS * t17 + SIN ) * x( 0 ) * t23 * PW * t17;
  }

  bool check_x( Vector const & x ) const { return x( 1 ) > 0; }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    if ( check_x( x ) )
    {
      Vector g( 5 );
      f.setZero();
      for ( integer k = 0; k < 24; ++k )
      {
        real_type ff = f_val( x, k );
        f_grad( x, k, g );
        f += ff * g;
      }
    }
    else
    {
      f.fill( nan( "DeVilliersGlasser02" ) );
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    Vector g( 5 );
    Matrix h( 5, 5 ), J_full( 5, 5 );
    J_full.setZero();
    for ( integer k = 0; k < 24; ++k )
    {
      real_type ff = f_val( x, k );
      f_grad( x, k, g );
      f_hess( x, k, h );
      for ( integer ii = 0; ii < 5; ++ii )
        for ( integer jj = 0; jj < 5; ++jj ) J_full( ii, jj ) += ff * h( ii, jj ) + g( ii ) * g( jj );
    }
    J = J_full.sparseView();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 45, 2, 2.5, 1.5, 0.9;
  }

  virtual void check_if_admissible( Vector const & x ) const override
  {
    UTILS_ASSERT( x( 1 ) > 0, "DeVilliersGlasser02, x(1) = {} <= 0", x( 1 ) );
    UTILS_ASSERT( -500 <= x( 0 ) && x( 0 ) <= 500, "DeVilliersGlasser#2, x(0) = {} must be in [-500,500]", x( 0 ) );
    UTILS_ASSERT( -500 <= x( 1 ) && x( 1 ) <= 500, "DeVilliersGlasser#2, x(1) = {} must be in [-500,500]", x( 1 ) );
    UTILS_ASSERT( -500 <= x( 2 ) && x( 2 ) <= 500, "DeVilliersGlasser#2, x(2) = {} must be in [-500,500]", x( 2 ) );
    UTILS_ASSERT( -500 <= x( 3 ) && x( 3 ) <= 500, "DeVilliersGlasser#2, x(3) = {} must be in [-500,500]", x( 3 ) );
    UTILS_ASSERT( -500 <= x( 4 ) && x( 4 ) <= 500, "DeVilliersGlasser#2, x(4) = {} must be in [-500,500]", x( 4 ) );
  }

  virtual void bounding_box( Vector & L, Vector & U ) const override
  {
    L[0] = -500;
    U[0] = 500;
    L[1] = 0;
    U[1] = 500;
    L[2] = -500;
    U[2] = 500;
    L[3] = -500;
    U[3] = 500;
    L[4] = -500;
    U[4] = 500;
  }
};
