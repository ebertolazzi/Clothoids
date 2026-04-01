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

#define EXPONENTIAL_FUNCTION_BIBTEX                                        \
  "@article{LaCruz:2003,\n"                                                \
  "  author    = { William {La Cruz}  and  Marcos Raydan},\n"              \
  "  title     = {Nonmonotone Spectral Methods for Large-Scale Nonlinear " \
  "Systems},\n"                                                            \
  "  journal   = {Optimization Methods and Software},\n"                   \
  "  year      = {2003},\n"                                                \
  "  volume    = {18},\n"                                                  \
  "  number    = {5},\n"                                                   \
  "  pages     = {583--599},\n"                                            \
  "  publisher = {Taylor & Francis},\n"                                    \
  "  doi       = {10.1080/10556780310001610493},\n"                        \
  "}\n"

class ExponentialFunction1 : public NonlinearSystem
{
public:
  ExponentialFunction1( integer neq ) : NonlinearSystem( "Exponential Function N.1", EXPONENTIAL_FUNCTION_BIBTEX, neq )
  {
    check_min_equations( n, 1 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = exp( x( 0 ) - 1 ) - 1;
    for ( integer i = 1; i < n; ++i ) f( i ) = ( i + 1 ) * ( exp( x( i ) - 1 ) - x( i ) );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    J.insert( 0, 0 ) = exp( x( 0 ) - 1 );
    for ( integer i = 1; i < n; ++i ) J.insert( i, i ) = ( i + 1 ) * ( exp( x( i ) - 1 ) - 1 );
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( n / ( n - 1.0 ) );
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class ExponentialFunction2 : public NonlinearSystem
{
public:
  ExponentialFunction2( integer neq ) : NonlinearSystem( "Exponential Function N.2", EXPONENTIAL_FUNCTION_BIBTEX, neq )
  {
    check_min_equations( n, 1 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = exp( x( 0 ) ) - 1;
    for ( integer i = 1; i < n; ++i ) f( i ) = ( ( i + 1 ) / 10.0 ) * ( exp( x( i ) ) + x( i - 1 ) - 1 );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    J.insert( 0, 0 ) = exp( x( 0 ) );
    for ( integer i = 1; i < n; ++i )
    {
      J.insert( i, i )     = ( ( i + 1 ) / 10.0 ) * exp( x( i ) );
      J.insert( i, i - 1 ) = ( ( i + 1 ) / 10.0 );
    }
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 1.0 / ( n * n ) );
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class ExponentialFunction3 : public NonlinearSystem
{
public:
  ExponentialFunction3( integer neq ) : NonlinearSystem( "Exponential Function N.3", EXPONENTIAL_FUNCTION_BIBTEX, neq )
  {
    check_min_equations( n, 1 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( n - 1 ) = ( 0.1 * n ) * ( 1 - exp( -x( n - 1 ) * x( n - 1 ) ) );
    for ( integer i = 0; i < n - 1; ++i )
      f( i ) = ( 0.1 * ( i + 1 ) ) * ( 1 - x( i ) * x( i ) - exp( -x( i ) * x( i ) ) );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer i = 0; i < n - 1; ++i ) J.insert( i, i ) = 0.2 * ( i + 1 ) * x( i ) * ( exp( -x( i ) * x( i ) ) - 1 );
    J.insert( n - 1, n - 1 ) = 0.2 * n * x( n - 1 ) * exp( -x( n - 1 ) * x( n - 1 ) );
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    real_type bf = 1.0 / ( 4.0 * n * n );
    for ( integer i = 0; i < n; ++i ) x0( i ) = ( i + 1 ) * bf;
  }
};
