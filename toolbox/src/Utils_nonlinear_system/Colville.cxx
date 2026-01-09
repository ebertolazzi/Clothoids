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

class Colville : public NonlinearSystem
{
public:
  Colville()
    : NonlinearSystem(
        "Colville Polynomial",
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

    real_type x1_2 = x1 * x1;
    real_type x1_3 = x1_2 * x1;
    real_type x3_2 = x3 * x3;
    real_type x3_3 = x3_2 * x3;

    f( 0 ) = 400.0 * x1_3 - 400.0 * x2 * x1 + 2.0 * x1 - 2.0;
    f( 1 ) = -200.0 * x1_2 + 220.2 * x2 + 19.8 * x4 - 40.0;
    f( 2 ) = -360.0 * x3 * x4 + 360.0 * x3_3 + 2.0 * x3 - 2.0;
    f( 3 ) = 180.0 * x4 - 180.0 * x3_2 + 20.2 * x4 + 19.8 * x2 - 40.0;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    real_type x1 = x( 0 );
    real_type x2 = x( 1 );
    real_type x3 = x( 2 );
    real_type x4 = x( 3 );

    J.insert( 0, 0 ) = 1200.0 * x1 * x1 - 400.0 * x2 + 2.0;
    J.insert( 0, 1 ) = -400.0 * x1;

    J.insert( 1, 0 ) = -400.0 * x1;
    J.insert( 1, 1 ) = 220.2;
    J.insert( 1, 3 ) = 19.8;

    J.insert( 2, 2 ) = -360.0 * x4 + 1080.0 * x3 * x3 + 2.0;
    J.insert( 2, 3 ) = -360.0 * x3;

    J.insert( 3, 1 ) = 19.8;
    J.insert( 3, 2 ) = -360.0 * x3;
    J.insert( 3, 3 ) = 200.2;

    J.makeCompressed();
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
    x0 << 0.5, 1.0, -0.5, -1.0;
  }
};
