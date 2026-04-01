/*\
 |
 |  Author:
 |    Enrico Bertolazzi
 |    University of Trento
 |    Department of Industrial Engineering
 |    Via Sommarive 9, I-38123, Povo, Trento, Italy
 |    email: enrico.bertolazzi@unitn.it
\*/

// AZZZ QUALCOSA NON TORNA @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class SpedicatoFunction17 : public NonlinearSystem
{
public:
  SpedicatoFunction17( integer neq )
    : NonlinearSystem(
        "Spedicato N.17",
        "@Article{Spedicato1997,\n"
        "  author  = {Spedicato, E. and Huang, Z.},\n"
        "  title   = {Numerical experience with newton-like methods\n"
        "             for nonlinear algebraic systems},\n"
        "  journal = {Computing},\n"
        "  year    = {1997},\n"
        "  volume  = {58},\n"
        "  number  = {1},\n"
        "  pages   = {69--89},\n"
        "  doi     = {10.1007/BF02684472},\n"
        "}\n\n"
        "@book{meresoo:1990,\n"
        "  title     = {Test Examples of Systems of Nonlinear Equations: "
        "Version 3-90},\n"
        "  author    = {Meresoo, T. and Roose, A. and Kulla,\n"
        "               V. and Estonian Software and Computer Service "
        "Company},\n"
        "  year      = 1990,\n"
        "  publisher = {Estonian Software and Computer Service Company}\n"
        "}\n",
        neq )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 )     = x( 0 );
    f( n - 1 ) = x( n - 1 ) - 20;
    for ( integer i = 1; i < n - 1; ++i )
      f( i ) = x( i + 1 ) + x( i ) + x( i - 1 ) + power2( x( i + 1 ) - x( i - 1 ) ) / 4;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    J.insert( 0, 0 )         = 1;
    J.insert( n - 1, n - 1 ) = 1;
    for ( integer i = 1; i < n - 1; ++i )
    {
      real_type bf         = x( i + 1 ) - x( i - 1 );
      J.insert( i, i )     = 1;
      J.insert( i, i + 1 ) = 1 + bf / 2;
      J.insert( i, i - 1 ) = 1 - bf / 2;
    }
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 10 );
  }

  virtual void check_if_admissible( Vector const & x ) const override
  {
    for ( integer i = 0; i < n; ++i ) UTILS_ASSERT( std::abs( x( i ) ) < 10000, "Bad range" );
  }

  virtual void bounding_box( Vector & L, Vector & U ) const override
  {
    U.fill( 10000 );
    L.fill( -10000 );
  }
};
