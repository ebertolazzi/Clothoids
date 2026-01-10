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

#define SoniaKrzyworzcka1_BIBTEX                                      \
  "@article{Krzyworzcka:1996,\n"                                      \
  "  author  = {Sonia Krzyworzcka},\n"                                \
  "  title   = {Extension of the Lanczos and {CGS}\n"                 \
  "             methods to systems of nonlinear equations},\n"        \
  "  journal = {Journal of Computational and Applied Mathematics},\n" \
  "  volume  = {69},\n"                                               \
  "  number  = {1},\n"                                                \
  "  pages   = {181--190},\n"                                         \
  "  year    = {1996},\n"                                             \
  "  doi     = {10.1016/0377-0427(95)00032-1}\n"                      \
  "}\n"

class SoniaKrzyworzcka1 : public NonlinearSystem
{
public:
  SoniaKrzyworzcka1() : NonlinearSystem( "Sonia Krzyworzcka example 2", SoniaKrzyworzcka1_BIBTEX, 6 ) {}

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = -0.75 - 0.5 * x( 1 ) * x( 1 ) * x( 3 ) * x( 5 ) - x( 0 );
    f( 1 ) = -0.405 * exp( 1 + x( 0 ) * x( 1 ) ) + 1.405 - x( 1 );
    f( 2 ) = 0.5 * x( 3 ) * x( 5 ) - 1.5 - x( 2 );
    f( 3 ) = 0.605 * exp( 1 - x( 2 ) * x( 2 ) ) + 0.395 - x( 3 );
    f( 4 ) = 0.5 * x( 1 ) * x( 5 ) - 1.5 - x( 4 );
    f( 5 ) = x( 0 ) * x( 5 ) - x( 5 );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    J.insert( 0, 0 ) = -1;
    J.insert( 0, 1 ) = -1.0 * x( 1 ) * x( 3 ) * x( 5 );
    J.insert( 0, 3 ) = -0.5 * x( 1 ) * x( 1 ) * x( 5 );
    J.insert( 0, 5 ) = -0.5 * x( 1 ) * x( 1 ) * x( 3 );

    J.insert( 1, 0 ) = -0.405 * x( 1 ) * exp( x( 0 ) * x( 1 ) + 1 );
    J.insert( 1, 1 ) = -0.405 * x( 0 ) * exp( x( 0 ) * x( 1 ) + 1 ) - 1;

    J.insert( 2, 2 ) = -1;
    J.insert( 2, 3 ) = 0.5 * x( 5 );
    J.insert( 2, 5 ) = 0.5 * x( 3 );

    J.insert( 3, 2 ) = -1.210 * x( 2 ) * exp( -x( 2 ) * x( 2 ) + 1 );
    J.insert( 3, 3 ) = -1;

    J.insert( 4, 1 ) = 0.5 * x( 5 );
    J.insert( 4, 4 ) = -1;
    J.insert( 4, 5 ) = 0.5 * x( 1 );

    J.insert( 5, 0 ) = x( 5 );
    J.insert( 5, 5 ) = x( 0 ) - 1;

    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << -1, 1, -1, 1, -1, 1;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 0.1, 0.1, 0.1, 0.1, 0.1, 0.1;
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class SoniaKrzyworzcka2 : public NonlinearSystem
{
public:
  SoniaKrzyworzcka2() : NonlinearSystem( "Sonia Krzyworzcka example 3", SoniaKrzyworzcka1_BIBTEX, 10 ) {}

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    for ( integer i = 0; i < n; ++i ) f( i ) = ( 3 - 5 * x( i ) ) * x( i ) + 1;
    for ( integer i = 0; i < n - 1; ++i ) f( i ) -= 2 * x( i + 1 );
    for ( integer i = 1; i < n; ++i ) f( i ) -= x( i - 1 );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer i = 0; i < n; ++i ) J.insert( i, i ) = 3 - 10 * x( i );
    for ( integer i = 0; i < n - 1; ++i ) J.insert( i, i + 1 ) = -2;
    for ( integer i = 1; i < n; ++i ) J.insert( i, i - 1 ) = -1;
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( -1 );
  }
};
