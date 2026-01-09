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

class BrownAlmostLinearFunction : public NonlinearSystem
{
  using Matrix = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;

public:
  BrownAlmostLinearFunction( integer neq )
    : NonlinearSystem(
        "Brown almost linear function",
        "@article{Brown:1968,\n"
        "  author    = {Brown, Kenneth M.},\n"
        "  title     = {A Quadratically Convergent Newton-Like Method "
        "Based\n"
        "               Upon Gaussian-Elimination},\n"
        "  volume    = {6},\n"
        "  number    = {4},\n"
        "  year      = {1969},\n"
        "  pages     = {560--569}\n"
        "  publisher = {SIAM Journal on Numerical Analysis},\n"
        "}\n\n"
        "@article{More:1981,\n"
        "  author  = {Mor{\'e}, Jorge J. and Garbow, Burton S. and "
        "Hillstrom, Kenneth E.},\n"
        "  title   = {Testing Unconstrained Optimization Software},\n"
        "  journal = {ACM Trans. Math. Softw.},\n"
        "  volume  = {7},\n"
        "  number  = {1},\n"
        "  year    = {1981},\n"
        "  pages   = {17--41},\n"
        "  doi     = {10.1145/355934.355936},\n"
        "}\n",
        neq )
  {
    check_min_equations( neq, 2 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type sumx  = x.sum();
    real_type prodx = x.prod();
    for ( integer i = 0; i < n - 1; ++i ) f( i ) = x( i ) + ( sumx - ( n + 1 ) );
    f( n - 1 ) = prodx - 1;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    Matrix J_full( n, n );
    J_full.fill( 1 );
    for ( integer i = 0; i < n - 1; ++i ) J_full( i, i ) = 2;

    // last row
    for ( integer j = 0; j < n; ++j )
    {
      real_type prod = 1;
      for ( integer k = 0; k < n; ++k )
      {
        if ( k != j ) prod *= x( k );
      }
      J_full( n - 1, j ) = prod;
    }
    J.resize( n, n );
    J = J_full.sparseView();
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
    x0.fill( 0.5 );
  }
};
