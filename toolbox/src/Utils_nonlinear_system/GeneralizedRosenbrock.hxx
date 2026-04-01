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

class GeneralizedRosenbrock : public NonlinearSystem
{
  real_type N;

public:
  GeneralizedRosenbrock( integer neq )
    : NonlinearSystem(
        "Generalized Rosenbrock function",
        "@article{Rosenbrock:1960,\n"
        "  author  = {Rosenbrock, H. H.},\n"
        "  title   = {An Automatic Method for Finding the Greatest\n"
        "             or Least Value of a Function},\n"
        "  journal = {The Computer Journal},\n"
        "  year    = {1960},\n"
        "  volume  = {3},\n"
        "  number  = {3},\n"
        "  pages   = {175--184},\n"
        "  doi = {10.1093/comjnl/3.3.175},\n"
        "}\n\n"
        "@article{More:1981,\n"
        "  author = {Mor{\'e}, Jorge J. and Garbow, Burton S. and "
        "Hillstrom, Kenneth E.},\n"
        "  title = {Testing Unconstrained Optimization Software},\n"
        "  journal = {ACM Trans. Math. Softw.},\n"
        "  volume = {7},\n"
        "  number = {1},\n"
        "  month = mar,\n"
        "  year = {1981},\n"
        "  pages = {17--41},\n"
        "  doi = {10.1145/355934.355936},\n"
        "}\n",
        neq )
    , N( 100 )
  {
    check_even( n, 2 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = -4 * N * ( x( 1 ) - power2( x( 0 ) ) ) * x( 0 ) + 2 * x( 0 ) - 2;
    for ( integer i = 1; i < n - 1; ++i )
      f( i ) = 2 * N * ( x( i ) - power2( x( i - 1 ) ) ) - 4 * N * ( x( i + 1 ) - power2( x( i ) ) ) * x( i ) +
               2 * x( i ) - 2;
    f( n - 1 ) = 2 * N * ( x( n - 1 ) - power2( x( n - 2 ) ) );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    J.insert( 0, 0 ) = 8 * N * power2( x( 0 ) ) - 4 * N * ( x( 1 ) - power2( x( 0 ) ) ) + 2;
    J.insert( 0, 1 ) = -4 * N * x( 0 );
    for ( integer i = 1; i < n - 1; ++i )
    {
      J.insert( i, i - 1 ) = -4 * N * x( i - 1 );
      J.insert( i, i )     = 2 + ( 12 * x( i ) * x( i ) - 4 * x( i + 1 ) + 2 ) * N;
      J.insert( i, i + 1 ) = -4 * N * x( i );
    }
    J.insert( n - 1, n - 2 ) = -4 * N * x( n - 2 );
    J.insert( n - 1, n - 1 ) = 2 * N;
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
    for ( integer i{ 0 }; i < n; i += 2 )
    {
      x0( i )     = -1.2;
      x0( i + 1 ) = 1;
    }
  }
};
