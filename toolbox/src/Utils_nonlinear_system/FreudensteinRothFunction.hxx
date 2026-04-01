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

class FreudensteinRothFunction : public NonlinearSystem
{
public:
  FreudensteinRothFunction()
    : NonlinearSystem(
        "Freudenstein-Roth function",
        "@article{Freudenstein:1963,\n"
        "  author  = {Freudenstein, Ferdinand and Roth, Bernhard},\n"
        "  title   = {Numerical Solution of Systems of Nonlinear "
        "Equations},\n"
        "  journal = {J. ACM},\n"
        "  year    = {1963},\n"
        "  volume  = {10},\n"
        "  number  = {4},\n"
        "  pages   = {550--556},\n"
        "  doi     = {10.1145/321186.321200}\n"
        "}\n\n"
        "@article{More:1981,\n"
        "  author  = {Mor{\'e}, Jorge J. and Garbow, Burton S. and "
        "Hillstrom, Kenneth E.},\n"
        "  title   = {Testing Unconstrained Optimization Software},\n"
        "  journal = {ACM Trans. Math. Softw.},\n"
        "  year    = {1981},\n"
        "  volume  = {7},\n"
        "  number  = {1},\n"
        "  pages   = {17--41},\n"
        "  doi     = {10.1145/355934.355936},\n"
        "}\n",
        2 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = x( 0 ) - power3( x( 1 ) ) + 5 * power2( x( 1 ) ) - 2 * x( 1 ) - 13;
    f( 1 ) = x( 0 ) + power3( x( 1 ) ) + power2( x( 1 ) ) - 14 * x( 1 ) - 29;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    J.insert( 0, 0 ) = 1;
    J.insert( 0, 1 ) = -3 * power2( x( 1 ) ) + 10 * x( 1 ) - 2;
    J.insert( 1, 0 ) = 1;
    J.insert( 1, 1 ) = 3 * power2( x( 1 ) ) + 2 * x( 1 ) - 14;
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 5, 4;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 0.5, -2;
  }
};
