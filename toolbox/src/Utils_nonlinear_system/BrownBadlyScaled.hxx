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

class BrownBadlyScaled : public NonlinearSystem
{
public:
  BrownBadlyScaled()
    : NonlinearSystem(
        "Brown Badly Scaled Function",
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
    real_type x1 = x( 0 );
    real_type x2 = x( 1 );

    f( 0 ) = ( 2 * x1 * ( 1 + x2 * x2 ) - 4 * x2 ) - 2000000.0;
    f( 1 ) = ( 2 * x2 * ( 1 + x1 * x1 ) - 4 * x1 ) - 0.000004;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type x1 = x( 0 );
    real_type x2 = x( 1 );

    J.insert( 0, 0 ) = 2 * ( 1 + x2 * x2 );
    J.insert( 0, 1 ) = J.insert( 1, 0 ) = 4 * ( x1 * x2 - 1 );
    J.insert( 1, 1 )                    = 2 * ( 1 + x1 * x1 );
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 1e6, 2e-6;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 1, 1;
  }
};
