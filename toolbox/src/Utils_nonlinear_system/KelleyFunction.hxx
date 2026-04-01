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

class KelleyFunction : public NonlinearSystem
{
public:
  KelleyFunction()
    : NonlinearSystem(
        "Kelley Function",
        "@book{Dennis:1996,\n"
        "  author    = {Dennis, J. and Schnabel, R.},\n"
        "  title     = {Numerical Methods for Unconstrained\n"
        "               Optimization and Nonlinear Equations},\n"
        "  publisher = {Society for Industrial and Applied Mathematics},\n"
        "  year      = {1996},\n"
        "  doi       = {10.1137/1.9781611971200},\n"
        "}\n",
        2 )
  {
  }

  virtual void evaluate( Vector const & x_in, Vector & f ) const override
  {
    real_type x = x_in[0];
    real_type y = x_in[1];
    f( 0 )      = x * x + y * y - 2;
    f( 1 )      = exp( x - 1 ) + y * y - 2;
  }

  virtual void jacobian( Vector const & x_in, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type x      = x_in[0];
    real_type y      = x_in[1];
    J.insert( 0, 0 ) = 2 * x;
    J.insert( 0, 1 ) = 2 * y;
    J.insert( 1, 0 ) = exp( x - 1 );
    J.insert( 1, 1 ) = 2 * y;
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
    x0 << 2, 1e-6;
  }
};
