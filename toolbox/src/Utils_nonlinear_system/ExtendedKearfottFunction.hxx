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

class ExtendedKearfottFunction : public NonlinearSystem
{
public:
  ExtendedKearfottFunction()
    : NonlinearSystem(
        "Extended Kearfott Function",
        "@phdthesis{Kearfott:1977,\n"
        "  author    = {Kearfott, Ralph Baker},\n"
        "  title     = {Computing the Degree of Maps and a Generalized\n"
        "               Method of Bisection},\n"
        "  year      = {1977},\n"
        "  note      = {AAI7723103},\n"
        "  publisher = {The University of Utah},\n"
        "}\n\n"
        "@Article{Kearfott:1979,\n"
        "  author  = {Kearfott, Baker},\n"
        "  title   = {An efficient degree-computation method for a "
        "generalized\n"
        "             method of bisection},\n"
        "  journal = {Numerische Mathematik},\n"
        "  year    = {1979},\n"
        "  volume  = {32},\n"
        "  number  = {2},\n"
        "  pages   = {109--127},\n"
        "  doi     = {10.1007/BF01404868}\n"
        "}\n",
        7 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    for ( integer i = 0; i < n - 1; ++i ) f( i ) = x( i ) * x( i ) - x( i + 1 );
    f( n - 1 ) = x( n - 1 ) * x( n - 1 ) - x( 0 );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer i = 0; i < n - 1; ++i )
    {
      J.insert( i, i )     = 2 * x( i );
      J.insert( i, i + 1 ) = -1;
    }
    J.insert( n - 1, n - 1 ) = 2 * x( n - 1 );
    J.insert( n - 1, 0 )     = -1;
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 0.1 );
  }
};
