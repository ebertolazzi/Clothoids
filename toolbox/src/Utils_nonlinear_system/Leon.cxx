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

class Leon : public NonlinearSystem
{
public:
  Leon()
    : NonlinearSystem(
        "Leon cubic valley function",
        "@book{brent2013,\n"
        "  author    = {Brent, R.P.},\n"
        "  title     = {Algorithms for Minimization Without Derivatives},\n"
        "  isbn      = {9780486143682},\n"
        "  series    = {Dover Books on Mathematics},\n"
        "  year      = {2013},\n"
        "  publisher = {Dover Publications}\n"
        "}\n",
        2 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type x0_3 = x( 0 ) * x( 0 ) * x( 0 );
    f( 0 )         = -600.0 * ( x( 1 ) - x0_3 ) * x( 0 ) * x( 0 ) - 2.0 * ( 1.0 - x( 0 ) );
    f( 1 )         = 200.0 * ( x( 1 ) - x0_3 );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type x0_4   = x( 0 ) * x( 0 ) * x( 0 ) * x( 0 );
    J.insert( 0, 0 ) = -1200.0 * x( 0 ) * x( 1 ) + 3000.0 * x0_4 + 2.0;
    J.insert( 0, 1 ) = -600.0 * x( 0 ) * x( 0 );
    J.insert( 1, 0 ) = -600.0 * x( 0 ) * x( 0 );
    J.insert( 1, 1 ) = 200.0;
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
    x0 << -1.2, -1.0;
  }
};
