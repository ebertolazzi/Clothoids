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

class DeistSefor : public NonlinearSystem
{
  real_type beta[6];

public:
  // sum log(xi-2)^2+log(xi-10)^2 - prod( xi) ^(1/5)
  DeistSefor()
    : NonlinearSystem(
        "DeistSefor function",
        "@Article{Martinez1980,\n"
        "  author  = {Mart{\\'i}nez, Jos{\\'e} Mario},\n"
        "  title   = {Solving nonlinear simultaneous equations with\n"
        "             a generalization of Brent's method},\n"
        "  journal = {BIT Numerical Mathematics},\n"
        "  year    = {1980},\n"
        "  volume  = {20},\n"
        "  number  = {4},\n"
        "  pages   = {501--510},\n"
        "  doi     = {10.1007/BF01933643}\n"
        "}\n",
        6 )
  {
    beta[0] = 0.02249;
    beta[1] = 0.02166;
    beta[2] = 0.02083;
    beta[3] = 0.02;
    beta[4] = 0.01918;
    beta[5] = 0.01835;
  }

  real_type cot( real_type x ) const { return 1 / tan( x ); }
  real_type csc2( real_type x ) const { return 1 / power2( sin( x ) ); }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    for ( integer i = 0; i < n; ++i )
    {
      f( i ) = 0;
      for ( integer j = 0; j < n; ++j )
        if ( i != j ) f( i ) += cot( beta[i] * x( j ) );
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer i = 0; i < n; ++i )
    {
      for ( integer j = 0; j < n; ++j )
        if ( i != j ) J.insert( i, j ) = -beta[i] / power2( sin( beta[i] * x( j ) ) );
    }
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 75 );
  }
};
