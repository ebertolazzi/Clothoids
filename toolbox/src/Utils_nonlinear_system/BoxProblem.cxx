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

class BoxProblem : public NonlinearSystem
{
public:
  BoxProblem()
    : NonlinearSystem(
        "Box Problem",
        "@article{Box:1966,\n"
        "  author  = {Box, M. J.},\n"
        "  title   = {A Comparison of Several Current Optimization "
        "Methods,\n"
        "             and the use of Transformations in Constrained "
        "Problems},\n"
        "  journal = {The Computer Journal},\n"
        "  volume  = {9},\n"
        "  number  = {1},\n"
        "  pages   = {67-77},\n"
        "  year    = {1966},\n"
        "  doi     = {10.1093/comjnl/9.1.67},\n"
        "}\n\n"
        "@article{More:1981,\n"
        "  author  = {Mor{\'e}, Jorge J. and Garbow, Burton S. and "
        "Hillstrom, Kenneth E.},\n"
        "  title   = {Testing Unconstrained Optimization Software},\n"
        "  journal = {ACM Trans. Math. Softw.},\n"
        "  volume  = {7},\n"
        "  number  = {1},\n"
        "  year    = {1981},\n"
        "  pages   = {17--41},\n"
        "  doi     = {10.1145/355934.355936},\n"
        "}\n",
        3 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type const tmp0{ exp( -0.1 ) - exp( -1 ) };
    real_type const tmp1{ exp( -0.2 ) - exp( -2 ) };
    real_type const tmp2{ exp( -0.3 ) - exp( -3 ) };

    f( 0 ) = exp( -0.1 * x( 0 ) ) - exp( -0.1 * x( 1 ) ) - x( 2 ) * tmp0;
    f( 1 ) = exp( -0.2 * x( 0 ) ) - exp( -0.2 * x( 1 ) ) - x( 2 ) * tmp1;
    f( 2 ) = exp( -0.3 * x( 0 ) ) - exp( -0.3 * x( 1 ) ) - x( 2 ) * tmp2;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    real_type const tmp0{ exp( -0.1 ) - exp( -1 ) };
    real_type const tmp1{ exp( -0.2 ) - exp( -2 ) };
    real_type const tmp2{ exp( -0.3 ) - exp( -3 ) };

    J.resize( n, n );
    J.setZero();

    J.insert( 0, 0 ) = -0.1 * exp( -0.1 * x( 0 ) );
    J.insert( 0, 1 ) = 0.1 * exp( -0.1 * x( 1 ) );
    J.insert( 0, 2 ) = -tmp0;

    J.insert( 1, 0 ) = -0.2 * exp( -0.2 * x( 0 ) );
    J.insert( 1, 1 ) = 0.2 * exp( -0.2 * x( 1 ) );
    J.insert( 1, 2 ) = -tmp1;

    J.insert( 2, 0 ) = -0.3 * exp( -0.3 * x( 0 ) );
    J.insert( 2, 1 ) = 0.3 * exp( -0.3 * x( 1 ) );
    J.insert( 2, 2 ) = -tmp2;

    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 2 );
    auto & x0{ x_vec[0] };
    auto & x1{ x_vec[1] };
    x0.resize( n );
    x1.resize( n );
    x0 << 1, 10, 1;
    x1 << 10, 1, -1;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 0, 1, 0;
  }
};
