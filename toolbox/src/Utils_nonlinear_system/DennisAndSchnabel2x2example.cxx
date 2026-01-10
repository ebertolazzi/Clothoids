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

class DennisAndSchnabel2x2example : public NonlinearSystem
{
public:
  DennisAndSchnabel2x2example()
    : NonlinearSystem(
        "Dennis and Schnabel 2 by 2 example",
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

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = x( 0 ) + x( 1 ) - 3;
    f( 1 ) = power2( x( 0 ) ) + power2( x( 1 ) ) - 9;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    J.insert( 0, 0 ) = 1;
    J.insert( 0, 1 ) = 1;
    J.insert( 1, 0 ) = 2 * x( 0 );
    J.insert( 1, 1 ) = 2 * x( 1 );
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 0, 3;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 1, 5;
  }
};
