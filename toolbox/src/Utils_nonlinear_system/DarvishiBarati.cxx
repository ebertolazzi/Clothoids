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

class DarvishiBarati : public NonlinearSystem
{
public:
  DarvishiBarati()
    : NonlinearSystem(
        "DarvishiBarati",
        "@article{Darvishi:2007,\n"
        "  author  = {Darvishi, M.T. and Barati, A.},\n"
        "  title   = {Super cubic iterative methods to solve systems\n"
        "             of nonlinear equations},\n"
        "  journal = {Applied Mathematics and Computation},\n"
        "  volume  = {188},\n"
        "  number  = {2},\n"
        "  pages   = {1678--1685},\n"
        "  year    = {2007},\n"
        "  doi     = {10.1016/j.amc.2006.11.022}\n"
        "}\n",
        2 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type x1 = x( 0 );
    real_type x2 = x( 1 );
    f( 0 )       = exp( x1 + x2 ) + x1 * cos( x2 );
    f( 1 )       = x1 + x2 - 1;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type x1     = x( 0 );
    real_type x2     = x( 1 );
    J.insert( 0, 0 ) = exp( x1 + x2 ) + cos( x2 );
    J.insert( 0, 1 ) = exp( x1 + x2 ) - x1 * sin( x2 );
    J.insert( 1, 0 ) = 1;
    J.insert( 1, 1 ) = 1;
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << -4, 5;
  }
};
