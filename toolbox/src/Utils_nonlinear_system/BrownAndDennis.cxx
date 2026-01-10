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

class BrownAndDennis : public NonlinearSystem
{
  using Matrix = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;

public:
  BrownAndDennis()
    : NonlinearSystem(
        "Brown and Dennis Function",
        "@book{brent2013,\n"
        "  author    = {Brent, R.P.},\n"
        "  title     = {Algorithms for Minimization Without Derivatives},\n"
        "  isbn      = {9780486143682},\n"
        "  series    = {Dover Books on Mathematics},\n"
        "  year      = {2013},\n"
        "  publisher = {Dover Publications}\n"
        "}\n",
        4 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type x1 = x( 0 );
    real_type x2 = x( 1 );
    real_type x3 = x( 2 );
    real_type x4 = x( 3 );

    f( 0 ) = f( 1 ) = f( 2 ) = f( 3 ) = 0;

    for ( integer i = 0; i < 20; ++i )
    {
      real_type c = ( i + 1.0 ) / 5.0;

      real_type sinc = sin( c );
      real_type f1   = x1 + c * x2 - exp( c );
      real_type f2   = x3 + sinc * x4 - cos( c );

      real_type f11 = f1 * f1;
      real_type f22 = f2 * f2;

      real_type tmp = f11 + f22;
      real_type t1  = f1 * tmp;
      real_type t2  = f2 * tmp;

      f( 0 ) += t1;
      f( 1 ) += t1 * c;
      f( 2 ) += t2;
      f( 3 ) += t2 * sinc;
    }

    f( 0 ) *= 4;
    f( 1 ) *= 4;
    f( 2 ) *= 4;
    f( 3 ) *= 4;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    real_type x1 = x( 0 );
    real_type x2 = x( 1 );
    real_type x3 = x( 2 );
    real_type x4 = x( 3 );

    Matrix J_full( n, n );

    J_full.setZero();

    for ( integer i = 0; i < 20; ++i )
    {
      real_type c    = ( i + 1.0 ) / 5.0;
      real_type sinc = sin( c );

      real_type f1 = x1 + c * x2 - exp( c );
      real_type f2 = x3 + sinc * x4 - cos( c );

      real_type f11 = f1 * f1;
      real_type f12 = f1 * f2;
      real_type f22 = f2 * f2;

      // f(0) += f1^3+f1*f2^2;
      real_type t1 = 3 * f11 + f22;
      real_type t2 = 2 * f12;
      J_full( 0, 0 ) += t1;
      J_full( 0, 1 ) += t1 * c;
      J_full( 0, 2 ) += t2;
      J_full( 0, 3 ) += t2 * sinc;

      // f(1) += (f1^3+f1*f2^2) * c;
      t1 *= c;
      t2 *= c;
      J_full( 1, 0 ) += t1;
      J_full( 1, 1 ) += t1 * c;
      J_full( 1, 2 ) += t2;
      J_full( 1, 3 ) += t2 * sinc;

      // f(2) += f1^2*f2 + f2^3;
      t1 = 2 * f12;
      t2 = f11 + 3 * f22;
      J_full( 2, 0 ) += t1;
      J_full( 2, 1 ) += t1 * c;
      J_full( 2, 2 ) += t2;
      J_full( 2, 3 ) += t2 * sinc;

      // f(3) += (f1^2*f2 + f2^3) * sinc;
      t1 *= sinc;
      t2 *= sinc;
      J_full( 3, 0 ) += t1;
      J_full( 3, 1 ) += t1 * c;
      J_full( 3, 2 ) += t2;
      J_full( 3, 3 ) += t2 * sinc;
    }
    J_full *= 4;
    J.resize( n, n );
    J = J_full.sparseView();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << -11.59443990476216538261421860886232167047, 13.20363005120720382131746247754292259854,
      -0.4034394881768595196441977938137369661485, 0.2367787744557362991471155778176237456915;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 25.0, 5.0, -5.0, -1.0;
  }
};
