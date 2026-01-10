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

class Powell3D : public NonlinearSystem
{
  using Matrix = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;

public:
  Powell3D()
    : NonlinearSystem(
        "Powell 3D Function",
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
  }

  virtual void evaluate( Vector const & X, Vector & f ) const override
  {
    real_type x = X( 0 );
    real_type y = X( 1 );
    real_type z = X( 2 );

    real_type t3  = x - y;
    real_type t6  = t3 * t3;
    real_type t8  = pow( 1.0 + t6, 2.0 );
    real_type t13 = 2.0 * t3 / t8;
    real_type t17 = cos( m_pi * y * z / 2.0 );
    real_type t18 = m_pi * t17;
    f( 0 )        = t13;
    f( 1 )        = -t13 - z * t18 / 2.0;
    f( 2 )        = -y * t18 / 2.0;

    if ( y != 0 )
    {
      real_type tt = exp( -( x + 2.0 * y + z ) / y ) / y;
      f( 0 ) -= tt;
      f( 1 ) += ( x + z ) * tt / y;
      f( 2 ) -= tt;
    }
  }

  virtual void jacobian( Vector const & X, SparseMatrix & J ) const override
  {
    real_type x = X( 0 );
    real_type y = X( 1 );
    real_type z = X( 2 );

    real_type t3  = x * x;
    real_type t7  = x * y;
    real_type t9  = y * y;
    real_type t13 = -6.0 * t3 + 12.0 * t7 - 6.0 * t9 + 2.0;
    real_type t15 = t3 - 2.0 * t7 + t9 + 1.0;
    real_type t16 = t15 * t15;
    real_type t18 = 1 / t16 / t15;
    real_type t20 = -t18 * t13;
    real_type t22 = pow( x - y, 2.0 );
    real_type t26 = pow( 1.0 + t22, 2.0 );
    real_type t29 = m_pi * m_pi;
    real_type t30 = z * z;
    real_type t34 = m_pi * y * z / 2.0;
    real_type t35 = sin( t34 );
    real_type t42 = cos( t34 );
    real_type t46 = ( y * t35 * m_pi * z - 2.0 * t42 ) * m_pi / 4.0;

    Matrix J_full( n, n );
    J_full.setZero();

    J_full( 0, 0 ) = t18 * t13;
    J_full( 0, 1 ) = t20;
    J_full( 0, 2 ) = 0.0;
    J_full( 1, 0 ) = t20;
    J_full( 1, 1 ) = -8.0 * t18 * t22 + 2.0 / t26 + t35 * t30 * t29 / 4.0;
    J_full( 1, 2 ) = t46;
    J_full( 2, 0 ) = 0.0;
    J_full( 2, 1 ) = t46;
    J_full( 2, 2 ) = t35 * t9 * t29 / 4.0;

    if ( y != 0 )
    {
      real_type t2  = 2.0 * y;
      real_type t8  = exp( ( -x - t2 - z ) / y );
      real_type t10 = y * y;
      real_type t12 = t8 / t10;
      real_type t17 = t8 * ( -y + x + z ) / t10 / y;
      real_type t21 = t10 * t10;
      J_full( 0, 0 ) += t12;
      J_full( 0, 1 ) += -t17;
      J_full( 0, 2 ) += t12;
      J_full( 1, 0 ) += -t17;
      J_full( 1, 1 ) += ( x + z ) * t8 * ( x + z - t2 ) / t21;
      J_full( 1, 2 ) += -t17;
      J_full( 2, 0 ) += t12;
      J_full( 2, 1 ) += -t17;
      J_full( 2, 2 ) += t12;
    }
    J.resize( n, n );
    J = J_full.sparseView();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 1 );
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 0, 1, 2;
  }
};
