/*\
 |
 |  Author:
 |    Enrico Bertolazzi
 |    University of Trento
 |    Department of Industrial Engineering
 |    Via Sommarive 9, I-38123, Povo, Trento, Italy
 |    email: enrico.bertolazzi@unitn.it
\*/

// TEST 221
/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class VariablyDimensionedFunction : public NonlinearSystem
{
  mutable Vector gf2;

public:
  VariablyDimensionedFunction( integer neq )
    : NonlinearSystem(
        "Variably dimensioned function",
        "@book{brent2013,\n"
        "  author    = {Brent, R.P.},\n"
        "  title     = {Algorithms for Minimization Without Derivatives},\n"
        "  isbn      = {9780486143682},\n"
        "  series    = {Dover Books on Mathematics},\n"
        "  year      = {2013},\n"
        "  publisher = {Dover Publications}\n"
        "}\n"
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
        neq )
  {
    check_min_equations( neq, 2 );
    gf2.resize( n );
  }

  //  f = f1 * f1 * ( 1.0 + f1 * f1 ) + f2;

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type sum1 = 0;
    for ( integer j = 0; j < n; ++j ) sum1 += ( j + 1 ) * ( x( j ) - 1 );
    for ( integer j = 0; j < n; ++j ) f( j ) = x( j ) - 1 + ( j + 1 ) * sum1 * ( 1 + 2 * power2( sum1 ) );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type sum1 = 0;
    for ( integer j = 0; j < n; ++j ) sum1 += ( j + 1 ) * ( x( j ) - 1 );

    for ( integer k = 0; k < n; ++k )
    {
      for ( integer j = 0; j < n; ++j )
      {
        real_type tmp = ( k + 1 ) * ( 1 + 6 * power2( sum1 ) ) * ( j + 1 );
        if ( j == k ) tmp += 1;
        J.insert( k, j ) = tmp;
      }
    }
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    for ( integer k{ 0 }; k < n; ++k ) x0( k ) = 1 - real_type( k + 1 ) / real_type( n );
  }
};
