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

class Weibull : public NonlinearSystem
{
  real_type     z[99], y[99];
  integer const NPT;
  using Matrix = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;

public:
  Weibull()
    : NonlinearSystem(
        "Weibull function",
        "@Article{Shan70,\n"
        "  Title   = {Conditioning of quasi-{N}ewton methods "
        "for function minimization},\n"
        "  Author  = {David F. Shanno},\n"
        "  Journal = {Mathematics of Computation},\n"
        "  Year    = {1970},\n"
        "  Number  = {111},\n"
        "  Pages   = {647--656},\n"
        "  Volume  = {24}\n"
        "}\n",
        3 )
    , NPT( 99 )
  {
    for ( integer i = 0; i < NPT; ++i )
    {
      z[i] = ( i + 1 ) * 0.01;
      y[i] = 25 + pow( 50 * log( 1 / z[i] ), 2. / 3. );
    }
  }

  void map( Vector const & x, Vector & eq ) const
  {
    for ( integer k = 0; k < NPT; ++k )
    {
      real_type y2 = y[k] - x( 2 );
      real_type g  = pow( y2, x( 1 ) ) / x( 0 );
      real_type f  = exp( -g );
      eq( k )      = f - z[k];
    }
  }

  void Grad_map( Vector const & x, integer k, Vector & G ) const
  {
    real_type y2 = y[k] - x( 2 );
    real_type g  = pow( y2, x( 1 ) ) / x( 0 );
    real_type fg = g * exp( -g );
    real_type lg = log( y2 );
    G( 0 )       = fg / x( 0 );
    G( 1 )       = -fg * lg;
    G( 2 )       = fg * x( 1 ) / y2;
  }

  void Hess_map( Vector const & x, integer k, Matrix & H ) const
  {
    real_type y2   = y[k] - x( 2 );
    real_type g    = pow( y2, x( 1 ) ) / x( 0 );
    real_type fg   = g * exp( -g );
    real_type fg_1 = ( 1 - g ) * exp( -g );
    real_type lg   = log( y2 );

    real_type g_x0 = -g / x( 0 );
    real_type g_x1 = g * lg;
    real_type g_x2 = -g * x( 1 ) / y2;

    // real_type lg_x2 = -1/y2;

    H( 0, 0 ) = ( fg_1 * g_x0 - fg / x( 0 ) ) / x( 0 );
    H( 0, 1 ) = fg_1 * g_x1 / x( 0 );
    H( 0, 2 ) = fg_1 * g_x2 / x( 0 );
    ;

    H( 1, 0 ) = H( 0, 1 );
    H( 1, 1 ) = -fg_1 * g_x1 * lg;
    H( 1, 2 ) = fg / y2 - fg_1 * g_x2 * lg;

    H( 2, 0 ) = H( 0, 2 );
    H( 2, 1 ) = H( 1, 2 );
    H( 2, 2 ) = ( fg_1 * g_x2 + fg / y2 ) * x( 1 ) / y2;
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

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 250, 0.3, 5;
  }
};
