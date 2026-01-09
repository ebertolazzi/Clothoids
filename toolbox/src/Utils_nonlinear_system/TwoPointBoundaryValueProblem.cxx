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

class TwoPointBoundaryValueProblem : public NonlinearSystem
{
public:
  TwoPointBoundaryValueProblem( integer neq )
    : NonlinearSystem(
        "Two-Point Boundary Value Problem",
        "@article{More:1979,\n"
        "  author  = {Mor{\'e}, Jorge J. and Cosnard, Michel Y.},\n"
        "  title   = {Numerical Solution of Nonlinear Equations},\n"
        "  journal = {ACM Trans. Math. Softw.},\n"
        "  year    = {1979},\n"
        "  volume  = {5},\n"
        "  number  = {1},\n"
        "  pages   = {64--85},\n"
        "  doi     = {10.1145/355815.355820},\n"
        "}\n\n",
        neq )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type h = 1.0 / ( n - 1.0 );

    f( 0 )     = x( 0 );
    f( n - 1 ) = x( n - 1 );
    for ( integer i = 1; i < n - 1; ++i )
    {
      real_type t = h * i;
      f( i )      = 2 * x( i ) - x( i - 1 ) - x( i + 1 ) + ( h * h / 2 ) * power3( x( i ) + t + 1 );
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    real_type h              = 1.0 / ( n - 1.0 );
    J.insert( 0, 0 )         = 1;
    J.insert( n - 1, n - 1 ) = 1;
    for ( integer i = 1; i < n - 1; ++i )
    {
      real_type t          = h * i;
      J.insert( i, i - 1 ) = -1;
      J.insert( i, i )     = 2 + 1.5 * ( h * h ) * power2( x( i ) + t + 1 );
      J.insert( i, i + 1 ) = -1;
    }
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    real_type h = 1.0 / ( n - 1.0 );
    for ( integer i{ 0 }; i < n; ++i )
    {
      real_type t = h * i;
      x0( i )     = t * ( t - 1 );
    }
  }
};
