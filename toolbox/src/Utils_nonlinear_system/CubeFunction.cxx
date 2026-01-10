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

class CubeFunction : public NonlinearSystem
{
public:
  CubeFunction()
    : NonlinearSystem(
        "Cube Function",
        "@inbook{Leon:1966,\n"
        "  title     = {Recent advances in optimization techniques: "
        "proceedings},\n"
        "  chapter   = {A comparison Among Eight Known Optimizing "
        "Procedures},\n"
        "  author    = {Leon, A.},\n"
        "  editor    = { Lavi, A. and Vogl, T.P.},\n"
        "  year      = {1966},\n"
        "  pages     = {28--46’,\n"
        "  publisher = {Wiley}\n"
        "}\n\n"
        "@article{doi:10.1137/0723046,\n"
        "  author  = {Grippo, L. and Lampariello, F. and Lucidi, S.},\n"
        "  title   = {A Nonmonotone Line Search Technique for Newton’s "
        "Method},\n"
        "  journal = {SIAM Journal on Numerical Analysis},\n"
        "  year    = {1986},\n"
        "  volume  = {23},\n"
        "  number  = {4},\n"
        "  pages   = {707--716},\n"
        "  doi     = {10.1137/0723046},\n"
        "}\n",
        2 )
  {
  }

  virtual void evaluate( Vector const & xx, Vector & f ) const override
  {
    real_type x  = xx( 0 );
    real_type y  = xx( 1 );
    real_type x3 = x * x * x;
    f( 0 )       = ( ( 600 * ( x3 - y ) * x ) + 2 ) * x - 2;
    f( 1 )       = 200 * ( y - x3 );
  }

  virtual void jacobian( Vector const & xx, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type x      = xx( 0 );
    real_type y      = xx( 1 );
    real_type x2     = x * x;
    real_type x3     = x2 * x;
    J.insert( 0, 0 ) = ( 3000 * x3 - 1200 * y ) * x + 2;
    J.insert( 0, 1 ) = J.insert( 1, 0 ) = -600 * x2;
    J.insert( 1, 1 )                    = 200;
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 1 );
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << -1.2, -1;
  }
};
