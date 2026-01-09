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

class WatsonFunction : public NonlinearSystem
{
  real_type t[29][31];
  using Matrix = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;

public:
  WatsonFunction()
    : NonlinearSystem(
        "Watson function",
        "@book{brent2013,\n"
        "  author    = {Brent, R.P.},\n"
        "  title     = {Algorithms for Minimization Without Derivatives},\n"
        "  isbn      = {9780486143682},\n"
        "  series    = {Dover Books on Mathematics},\n"
        "  year      = {2013},\n"
        "  publisher = {Dover Publications}\n"
        "}\n\n"
        "@book{kowalik1968methods,\n"
        "  author = {Kowalik, J.S. and Osborne, M.R.},\n"
        "  title  = {Methods for unconstrained optimization problems},\n"
        "  series = {Mathematical Linguistics and Automatic Language "
        "Processing},\n"
        "  year   = {1968},\n"
        "  publisher={American Elsevier Pub. Co.}\n"
        "}\n\n"
        "@article{More:1981,\n"
        "  author  = {Mor{\'e}, Jorge J. and Garbow, Burton S. and "
        "Hillstrom, Kenneth E.},\n"
        "  title   = {Testing Unconstrained Optimization Software},\n"
        "  journal = {ACM Trans. Math. Softw.},\n"
        "  year    = {1981},\n"
        "  volume  = {7},\n"
        "  number  = {1},\n"
        "  pages   = {17--41},\n"
        "  doi     = {10.1145/355934.355936},\n"
        "}\n",
        31 )
  {
    for ( integer i = 0; i < 29; ++i )
    {
      real_type ti = ( i + 1 ) / 29.0;
      t[i][0]      = 1;
      for ( integer j = 1; j < 31; ++j ) { t[i][j] = t[i][j - 1] * ti; }
    }
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f.setZero();
    for ( integer i = 0; i < 29; ++i )
    {
      real_type const * ti = t[i];
      real_type         fi = 0;
      for ( integer j = 1; j < 31; ++j ) fi += j * x( j ) * ti[j - 1];
      real_type tmp = 0;
      for ( integer j = 0; j < 31; ++j ) tmp += x( j ) * ti[j];
      f( i ) = fi - tmp * tmp - 1;
    }
    f( 29 ) = x( 0 );
    f( 30 ) = x( 1 ) - x( 0 ) * x( 0 ) - 1;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    Matrix J_full( n, n );
    J_full.setZero();
    for ( integer i = 0; i < 29; ++i )
    {
      real_type const * ti = t[i];
      for ( integer j = 1; j < 31; ++j ) J_full( i, j ) += j * ti[j - 1];
      real_type tmp = 0;
      for ( integer j = 0; j < 31; ++j ) tmp += x( j ) * ti[j];
      for ( integer j = 0; j < 31; ++j ) J_full( i, j ) -= 2 * tmp * ti[j];
    }
    J_full( 29, 0 ) = 1;
    J_full( 30, 0 ) = -2 * x( 0 );
    J_full( 30, 1 ) = 1;
    J.resize( n, n );
    J = J_full.sparseView();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.setZero();
  }
};
