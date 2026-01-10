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

class InfRefluxFunction : public NonlinearSystem
{
public:
  InfRefluxFunction()
    : NonlinearSystem(
        "InfReflux function",
        "@article{Paterson:1986,\n"
        "  author  = {W.R. Paterson},\n"
        "  title   = {A new method for solving a class of "
        "nonlinear equations},\n"
        "  journal = {Chemical Engineering Science},\n"
        "  year    = {1986},\n"
        "  volume  = {41},\n"
        "  number  = {7},\n"
        "  pages   = {1935--1937},\n"
        "  doi     = {10.1016/0009-2509(86)87077-4}\n"
        "}\n",
        1 )
  {
  }

  virtual void evaluate( Vector const & x_in, Vector & f ) const override
  {
    real_type x    = x_in[0];
    real_type arg1 = 1. / ( 1. - x );
    real_type arg2 = 0.95 - x;
    if ( x > 0 && arg1 > 0 && arg2 > 0 )
      f( 0 ) = ( 1. / 63. ) * log( x ) + ( 64. / 63. ) * log( arg1 ) + log( arg2 ) - log( 0.9 );
    else
      f( 0 ) = nan( "InfRefluxFunction" );
  }

  virtual void jacobian( Vector const & x_in, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    real_type x    = x_in[0];
    real_type arg1 = 1. / ( 1. - x );
    real_type arg2 = 0.95 - x;
    if ( x > 0 && arg1 > 0 && arg2 > 0 )
      J.insert( 0, 0 ) = ( 1. / 63. ) / x + ( 64. / 63. ) / ( 1. - x ) - 1. / ( 0.95 - x );
    else
      J.insert( 0, 0 ) = nan( "InfRefluxFunction" );
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 4 );
    auto & x0{ x_vec[0] };
    auto & x1{ x_vec[1] };
    auto & x2{ x_vec[2] };
    auto & x3{ x_vec[3] };
    x0.resize( n );
    x1.resize( n );
    x2.resize( n );
    x3.resize( n );

    x0 << 0.23;   // default initial guess (critical as close to J=0 for x=0.229)
    x1 << 0.228;  // an even harder initial point
    x2 << 0.6;    // and two relatively easy points
    x3 << 0.01;
  }

  virtual void check_if_admissible( Vector const & x_in ) const override
  {
    real_type x = x_in[0];
    UTILS_ASSERT( x > 0 && x < 0.95, "ARGUMENT ERROR" );
  }

  virtual void bounding_box( Vector & L, Vector & U ) const override
  {
    U[0] = 0.95;
    L[0] = 0;
  }
};
