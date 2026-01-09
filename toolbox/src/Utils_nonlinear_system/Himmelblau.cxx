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

class Himmelblau : public NonlinearSystem
{
public:
  Himmelblau()
    : NonlinearSystem(
        "Himmelblau function",
        "@book{himmelblau:1972,\n"
        "  author    = {Himmelblau, D.M.},\n"
        "  title     = {Applied nonlinear programming},\n"
        "  year      = {1972},\n"
        "  publisher = {McGraw-Hill}\n"
        "}\n\n"
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
    real_type x0_2 = x( 0 ) * x( 0 );
    real_type x1_2 = x( 1 ) * x( 1 );
    f( 0 )         = 4 * ( x0_2 + x( 1 ) - 11 ) * x( 0 ) + 2 * ( x( 0 ) + x1_2 - 7 );
    f( 1 )         = 2 * ( x0_2 + x( 1 ) - 11 ) + 4 * ( x( 0 ) + x1_2 - 7 ) * x( 1 );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    J.insert( 0, 0 ) = 8 * x( 0 ) * x( 0 ) + 4 * ( x( 0 ) * x( 0 ) + x( 1 ) - 11 ) + 2;
    J.insert( 0, 1 ) = 4 * x( 0 ) + 4 * x( 1 );
    J.insert( 1, 0 ) = 4 * x( 0 ) + 4 * x( 1 );
    J.insert( 1, 1 ) = 2 + 8 * x( 1 ) * x( 1 ) + 4 * ( x( 0 ) + x( 1 ) * x( 1 ) - 7 );
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 3.0, 2.0;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << -1.3, 2.7;
  }
};
