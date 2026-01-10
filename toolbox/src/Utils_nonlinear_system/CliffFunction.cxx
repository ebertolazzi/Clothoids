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

class CliffFunction : public NonlinearSystem
{
public:
  CliffFunction()
    : NonlinearSystem(
        "Cliff Function",
        "@article{Grippo:1991,\n"
        "  author  = {Grippo, L. and Lampariello, F. and Lucidi, S.},\n"
        "  title   = {A Class of Nonmonotone Stabilization Methods\n"
        "             in Unconstrained Optimization},\n"
        "  journal = {Numer. Math.},\n"
        "  year    = {1991},\n"
        "  volume  = {59},\n"
        "  number  = {1},\n"
        "  pages   = {779--805},\n"
        "  doi     = {10.1007/BF01385810},\n"
        "}\n",
        2 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type tmp = 1 - 20.0 * exp( 20.0 * ( x( 0 ) - x( 1 ) ) );
    f( 0 )        = x( 0 ) / 5000.0 - 3.0 / 5000.0 - tmp;
    f( 1 )        = tmp;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type tmp    = 400 * exp( 20.0 * ( x( 0 ) - x( 1 ) ) );
    J.insert( 0, 0 ) = 1 / 5000.0 + tmp;
    J.insert( 0, 1 ) = -tmp;
    J.insert( 1, 0 ) = -tmp;
    J.insert( 1, 1 ) = tmp;
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 3, 3 + log( 20.0 ) / 20.0;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 0, -1;
  }
};
