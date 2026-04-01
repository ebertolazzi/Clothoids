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

class DixonFunction : public NonlinearSystem
{
public:
  DixonFunction( integer neq )
    : NonlinearSystem(
        "Dixon Function",
        "@Article{Dixon1988,\n"
        "  author  = {Dixon, L. C. W. and Price, R. C.},\n"
        "  title   = {Numerical experience with the truncated Newton\n"
        "             method for unconstrained optimization},\n"
        "  journal = {Journal of Optimization Theory and Applications},\n"
        "  year    = {1988},\n"
        "  volume  = {56},\n"
        "  number  = {2},\n"
        "  pages   = {245--255},\n"
        "  doi     = {10.1007/BF00939410}\n"
        "}\n",
        neq )
  {
    check_min_equations( n, 2 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = 2 * ( x( 0 ) - 1 );
    for ( integer i = 1; i < n - 1; ++i )
      f( i ) = 8 * ( i + 1 ) * ( 2 * x( i ) * x( i ) - x( i - 1 ) ) -
               2 * ( i + 2 ) * ( 2 * x( i + 1 ) * x( i + 1 ) - x( i ) );
    f( n - 1 ) = 8 * n * ( 2 * x( n - 1 ) * x( n - 1 ) - x( n - 2 ) ) * x( n - 1 );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    J.insert( 0, 0 ) = 2;
    for ( integer i = 1; i < n - 1; ++i )
    {
      J.insert( i, i - 1 ) = -8 * ( i + 1 );
      J.insert( i, i )     = 32 * ( i + 1 ) * x( i ) + 2 * ( i + 2 );
      J.insert( i, i + 1 ) = -8 * ( i + 2 ) * x( i + 1 );
    }
    J.insert( n - 1, n - 2 ) = -8 * n * x( n - 1 );
    J.insert( n - 1, n - 1 ) = 32 * n * x( n - 1 ) * x( n - 1 ) + 8 * n * ( 2 * x( n - 1 ) * x( n - 1 ) - x( n - 2 ) );
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 1 );
  }
};
