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

class SIRtest : public NonlinearSystem
{
public:
  SIRtest( integer neq_in )
    : NonlinearSystem(
        "Semi-implicit approach Example 5",
        "@article{Scheffel:2009,\n"
        "  author  = {Jan Scheffel and Cristian Håkansson},\n"
        "  title   = {Solution of systems of nonlinear equations\n"
        "              – a semi-implicit approach},\n"
        "  journal = {Applied Numerical Mathematics},\n"
        "  volume  = {59},\n"
        "  number  = {10},\n"
        "  pages   = {2430--2443},\n"
        "  year    = {2009},\n"
        "  doi     = {10.1016/j.apnum.2009.05.002},\n"
        "}\n",
        neq_in )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    for ( integer i = 0; i < n - 1; ++i ) f( i ) = x( i ) - cos( x( i + 1 ) );
    f( n - 1 ) = x( n - 1 ) - 3 * cos( x( 0 ) );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer i = 0; i < n - 1; ++i )
    {
      J.insert( i, i )     = 1;
      J.insert( i, i + 1 ) = sin( x( i + 1 ) );
    }
    J.insert( n - 1, n - 1 ) = 1;
    J.insert( n - 1, 0 )     = 3 * sin( x( 0 ) );
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    if ( n == 2 )
      x0.fill( -2.0 );
    else
      x0.fill( 3.0 );
  }

  virtual void check_if_admissible( Vector const & x ) const override
  {
    for ( integer i = 0; i < n; ++i ) UTILS_ASSERT( x( i ) > -5 && x( i ) < 5, "Bad range" );
  }

  virtual void bounding_box( Vector & L, Vector & U ) const override
  {
    U.fill( 5 );
    L.fill( -5 );
  }
};
