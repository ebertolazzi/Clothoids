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

class Gauss : public NonlinearSystem
{
  real_type y[15];

  using Matrix = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;

public:
  Gauss()
    : NonlinearSystem(
        "Gaussian function",
        "@book{brent2013,\n"
        "  author    = {Brent, R.P.},\n"
        "  title     = {Algorithms for Minimization Without Derivatives},\n"
        "  isbn      = {9780486143682},\n"
        "  series    = {Dover Books on Mathematics},\n"
        "  year      = {2013},\n"
        "  publisher = {Dover Publications}\n"
        "}\n",
        3 )
  {
    y[0]  = 0.0009;
    y[1]  = 0.0044;
    y[2]  = 0.0175;
    y[3]  = 0.0540;
    y[4]  = 0.1295;
    y[5]  = 0.2420;
    y[6]  = 0.3521;
    y[7]  = 0.3989;
    y[8]  = 0.3521;
    y[9]  = 0.2420;
    y[10] = 0.1295;
    y[11] = 0.0540;
    y[12] = 0.0175;
    y[13] = 0.0044;
    y[14] = 0.0009;
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type x1 = x( 0 );
    real_type x2 = x( 1 );
    real_type x3 = x( 2 );
    f( 0 ) = f( 1 ) = f( 2 ) = 0;
    for ( integer i = 0; i < 15; ++i )
    {
      real_type d1  = 0.5 * i;
      real_type d2  = 3.5 - d1 - x3;
      real_type arg = -0.5 * x2 * d2 * d2;
      real_type t   = x1 * exp( arg ) - y[i];

      f( 0 ) += 2.0 * exp( arg ) * t;
      f( 1 ) -= x1 * exp( arg ) * t * d2 * d2;
      f( 2 ) += 2.0 * x1 * x2 * exp( arg ) * t * d2;
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    real_type x1 = x( 0 );
    real_type x2 = x( 1 );
    real_type x3 = x( 2 );

    Matrix J_full( n, n );
    J.resize( n, n );
    J.setZero();
    J_full.setZero();

    for ( integer i = 0; i < 15; ++i )
    {
      real_type d1  = 0.5 * i;
      real_type d2  = 3.5 - d1 - x3;
      real_type arg = 0.5 * x2 * d2 * d2;
      real_type r   = exp( -arg );
      real_type t   = x1 * r - y[i];
      real_type t1  = 2.0 * x1 * r - y[i];

      real_type d2_2 = d2 * d2;
      real_type d2_4 = d2_2 * d2_2;

      J_full( 1 - 1, 1 - 1 ) += r * r;
      J_full( 2 - 1, 2 - 1 ) += r * t1 * d2_4;
      J_full( 3 - 1, 3 - 1 ) += r * ( x2 * t1 * d2_2 - t );
      J_full( 2 - 1, 1 - 1 ) -= r * t1 * d2_2;
      J_full( 3 - 1, 1 - 1 ) += d2 * r * t1;
      J_full( 3 - 1, 2 - 1 ) += d2 * r * ( t - arg * t1 );
    }

    J_full( 1 - 1, 1 - 1 ) *= 2.0;
    J_full( 2 - 1, 2 - 1 ) *= 0.5 * x1;
    J_full( 3 - 1, 3 - 1 ) *= 2.0 * x1 * x2;
    J_full( 3 - 1, 1 - 1 ) *= 2.0 * x2;
    J_full( 3 - 1, 2 - 1 ) *= 2.0 * x1;

    J_full( 1 - 1, 2 - 1 ) = J_full( 2 - 1, 1 - 1 );
    J_full( 1 - 1, 3 - 1 ) = J_full( 3 - 1, 1 - 1 );
    J_full( 2 - 1, 3 - 1 ) = J_full( 3 - 1, 2 - 1 );

    J.resize( n, n );
    J = J_full.sparseView();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 0.4, 1, 0;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.setZero();
  }
};
