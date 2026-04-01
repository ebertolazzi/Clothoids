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

class ExtendedPowellSingularFunction : public NonlinearSystem
{
  real_type const sqrt5;
  real_type const sqrt10;

public:
  ExtendedPowellSingularFunction()
    : NonlinearSystem(
        "Extended Powell singular function",
        "@article{Powell:1962,\n"
        "  author  = {Powell, M. J. D.},\n"
        "  title   = {An Iterative Method for Finding Stationary\n"
        "             Values of a Function of Several Variables},\n"
        "  journal = {The Computer Journal},\n"
        "  year    = {1962},\n"
        "  volume  = {5},\n"
        "  number  = {2},\n"
        "  pages   = {147--151},\n"
        "  doi     = {10.1093/comjnl/5.2.147}\n"
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
        4 )
    , sqrt5( sqrt( 5.0 ) )
    , sqrt10( sqrt( 10.0 ) )
  {
    check_four( n, 4 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    for ( integer i = 0; i < n; i += 4 )
    {
      f( i + 0 ) = x( i + 0 ) + 10 * x( i + 1 );
      f( i + 1 ) = sqrt5 * ( x( i + 2 ) - x( i + 3 ) );
      f( i + 2 ) = power2( x( i + 1 ) - 2 * x( i + 2 ) );
      f( i + 3 ) = sqrt10 * power2( x( i + 0 ) - x( i + 3 ) );
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer i = 0; i < n; i += 4 )
    {
      J.insert( i + 0, i + 0 ) = 1;
      J.insert( i + 0, i + 1 ) = 10;
      J.insert( i + 0, i + 2 ) = 0;
      J.insert( i + 0, i + 3 ) = 0;

      J.insert( i + 1, i + 0 ) = 0;
      J.insert( i + 1, i + 1 ) = 0;
      J.insert( i + 1, i + 2 ) = sqrt5;
      J.insert( i + 1, i + 3 ) = -sqrt5;

      J.insert( i + 2, i + 0 ) = 0;
      J.insert( i + 2, i + 1 ) = 2 * ( x( i + 1 ) - 2 * x( i + 2 ) );
      J.insert( i + 2, i + 2 ) = -4 * ( x( i + 1 ) - 2 * x( i + 2 ) );
      J.insert( i + 2, i + 3 ) = 0;

      J.insert( i + 3, i + 0 ) = 2 * sqrt10 * ( x( i + 0 ) - x( i + 3 ) );
      J.insert( i + 3, i + 1 ) = 0;
      J.insert( i + 3, i + 2 ) = 0;
      J.insert( i + 3, i + 3 ) = -2 * sqrt10 * ( x( i + 0 ) - x( i + 3 ) );
    }
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.setZero();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    for ( integer i{ 0 }; i < n; i += 4 )
    {
      x0( i + 0 ) = 3;
      x0( i + 1 ) = -1;
      x0( i + 2 ) = 0;
      x0( i + 3 ) = 1;
    }
  }
};
