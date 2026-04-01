/*\
 |
 |  Author:
 |    Enrico Bertolazzi
 |    University of Trento
 |    Department of Industrial Engineering
 |    Via Sommarive 9, I-38123, Povo, Trento, Italy
 |    email: enrico.bertolazzi@unitn.it
\*/

#define STRICT_CONVEX_FUNCTION_BIBTEX                                      \
  "@article{Raydan:1997,\n"                                                \
  "  author  = {Raydan, M.},\n"                                            \
  "  title   = {The Barzilai and Borwein Gradient Method for\n"            \
  "             the Large Scale Unconstrained Minimization Problem},\n"    \
  "  journal = {SIAM Journal on Optimization},\n"                          \
  "  volume  = {7},\n"                                                     \
  "  number  = {1},\n"                                                     \
  "  pages   = {26-33},\n"                                                 \
  "  year    = {1997},\n"                                                  \
  "  doi     = {10.1137/S1052623494266365},\n"                             \
  "}\n\n"                                                                  \
  "@article{LaCruz:2003,\n"                                                \
  "  author    = {William {La Cruz}  and  Marcos Raydan},\n"               \
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

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class StrictlyConvexFunction1 : public NonlinearSystem
{
public:
  StrictlyConvexFunction1( integer neq )
    : NonlinearSystem( "Strictly Convex Function 1", STRICT_CONVEX_FUNCTION_BIBTEX, neq )
  {
    check_min_equations( n, 1 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    for ( integer i = 0; i < n; ++i ) f( i ) = exp( x( i ) ) - 1;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer i = 0; i < n; ++i ) J.insert( i, i ) = exp( x( i ) );
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    for ( integer i{ 0 }; i < n; ++i ) x0( i ) = ( i + 1.0 ) / n;
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class StrictlyConvexFunction2 : public NonlinearSystem
{
public:
  StrictlyConvexFunction2( integer neq )
    : NonlinearSystem( "Strictly Convex Function 2", STRICT_CONVEX_FUNCTION_BIBTEX, neq )
  {
    check_min_equations( n, 1 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    for ( integer i = 0; i < n; ++i ) f( i ) = ( ( i + 1.0 ) / 10.0 ) * ( exp( x( i ) ) - 1 );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer i = 0; i < n; ++i ) J.insert( i, i ) = ( ( i + 1.0 ) / 10.0 ) * exp( x( i ) );
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
