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

class Box3 : public NonlinearSystem
{
  using Matrix = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;

public:
  Box3()
    : NonlinearSystem(
        "Box3",
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

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = f( 1 ) = f( 2 ) = 0;
    for ( integer i = 0; i < 10; ++i )
    {
      real_type c  = -( i + 1 ) / 10.0;
      real_type fi = exp( c * x( 0 ) ) - exp( c * x( 1 ) ) - x( 2 ) * ( exp( c ) - exp( 10 * c ) );

      real_type dfidx1 = c * exp( c * x( 0 ) );
      real_type dfidx2 = -c * exp( c * x( 1 ) );
      real_type dfidx3 = -( exp( c ) - exp( 10 * c ) );

      f( 0 ) += 2.0 * fi * dfidx1;
      f( 1 ) += 2.0 * fi * dfidx2;
      f( 2 ) += 2.0 * fi * dfidx3;
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    Matrix J_full( n, n );
    J_full.setZero();

    for ( integer i = 0; i < 10; ++i )
    {
      real_type c = -( i + 1 ) / 10.0;

      real_type fi = exp( c * x( 0 ) ) - exp( c * x( 1 ) ) - x( 2 ) * ( exp( c ) - exp( 10 * c ) );

      real_type dfidx1   = c * exp( c * x( 0 ) );
      real_type d2fidx11 = c * c * exp( c * x( 0 ) );
      real_type dfidx2   = -c * exp( c * x( 1 ) );
      real_type d2fidx22 = -c * c * exp( c * x( 1 ) );
      real_type dfidx3   = -( exp( c ) - exp( 10 * c ) );

      J_full( 0, 0 ) += 2.0 * dfidx1 * dfidx1 + 2.0 * fi * d2fidx11;
      J_full( 0, 1 ) += 2.0 * dfidx1 * dfidx2;
      J_full( 0, 2 ) += 2.0 * dfidx1 * dfidx3;

      J_full( 1, 0 ) += 2.0 * dfidx2 * dfidx1;
      J_full( 1, 1 ) += 2.0 * dfidx2 * dfidx2 + 2.0 * fi * d2fidx22;
      J_full( 1, 2 ) += 2.0 * dfidx2 * dfidx3;

      J_full( 2, 0 ) += 2.0 * dfidx3 * dfidx1;
      J_full( 2, 1 ) += 2.0 * dfidx3 * dfidx2;
      J_full( 2, 2 ) += 2.0 * dfidx3 * dfidx3;
    }
    J.resize( n, n );
    J = J_full.sparseView();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 1, 10, 1;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 0, 10, 5;
  }
};
