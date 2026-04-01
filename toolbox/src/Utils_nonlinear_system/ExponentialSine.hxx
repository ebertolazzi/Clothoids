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

class ExponentialSine : public NonlinearSystem
{
public:
  ExponentialSine()
    : NonlinearSystem(
        "Exponential sine",
        "@techreport{Nowak1991,\n"
        "  author = {U. Nowak and L. Weimann},\n"
        "  title  = {A Family of Newton Codes for Systems of "
        "Highly Nonlinear Equations},\n"
        "  number = {Technical Report TR-91-10},\n"
        "  year   = {1991}\n"
        "}\n",
        2 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type arg = x( 0 ) * x( 0 ) + x( 1 ) * x( 1 );
    f( 0 )        = exp( arg ) - 3.0;
    f( 1 )        = x( 0 ) + x( 1 ) - sin( 2.0 * ( x( 0 ) + x( 1 ) ) );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type arg    = x( 0 ) * x( 0 ) + x( 1 ) * x( 1 );
    real_type arg_0  = 2 * x( 0 );
    real_type arg_1  = 2 * x( 1 );
    J.insert( 0, 0 ) = exp( arg ) * arg_0;
    J.insert( 0, 1 ) = exp( arg ) * arg_1;
    J.insert( 1, 0 ) = J.insert( 1, 1 ) = 1 - 2 * cos( 2 * ( x( 0 ) + x( 1 ) ) );
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 5 );
    auto & x0{ x_vec[0] };
    auto & x1{ x_vec[1] };
    auto & x2{ x_vec[2] };
    auto & x3{ x_vec[3] };
    auto & x4{ x_vec[4] };
    x0.resize( n );
    x1.resize( n );
    x2.resize( n );
    x3.resize( n );
    x4.resize( n );
    x0 << 0.81, 0.82;
    x1 << 2.99714682530400, -2.95330710241260;
    x2 << 2.44169014107600, -2.51379895001340;
    x3 << 2.51913768310200, -2.85296012826840;
    x4 << -2.26027327648, -2.41295786268;
  }
};
