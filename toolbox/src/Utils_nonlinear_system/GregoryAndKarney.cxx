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

class GregoryAndKarney : public NonlinearSystem
{
public:
  GregoryAndKarney( integer n )
    : NonlinearSystem(
        "Gregory and Karney Tridiagonal Matrix Function",
        "@book{brent2013,\n"
        "  author    = {Brent, R.P.},\n"
        "  title     = {Algorithms for Minimization Without Derivatives},\n"
        "  isbn      = {9780486143682},\n"
        "  series    = {Dover Books on Mathematics},\n"
        "  year      = {2013},\n"
        "  publisher = {Dover Publications}\n"
        "}\n",
        n )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = x( 0 ) - x( 1 ) - 2;
    for ( integer i = 1; i < n - 1; ++i ) f( i ) = 2 * x( i ) - x( i - 1 ) - x( i + 1 );
    f( n - 1 ) = 2 * x( n - 1 ) - x( n - 2 );
  }

  virtual void jacobian( Vector const &, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer i = 0; i < n; ++i )
    {
      if ( i == 0 )
      {
        J.insert( i, i )     = 1;
        J.insert( i, i + 1 ) = -1;
      }
      else if ( i == n - 1 )
      {
        J.insert( i, i )     = 2;
        J.insert( i, i - 1 ) = -1;
      }
      else
      {
        J.insert( i, i - 1 ) = -1;
        J.insert( i, i )     = 2;
        J.insert( i, i + 1 ) = -1;
      }
    }
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    for ( integer i{ 0 }; i < n; ++i ) x0( i ) = 2.0 * ( n - i );
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.setZero();
  }
};
