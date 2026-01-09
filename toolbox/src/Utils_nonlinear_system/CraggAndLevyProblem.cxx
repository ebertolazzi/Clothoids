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

class CraggAndLevyProblem : public NonlinearSystem
{
public:
  CraggAndLevyProblem()
    : NonlinearSystem(
        "Cragg and Levy Problem",
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
        4 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = power2( exp( x( 0 ) ) - x( 1 ) );
    f( 1 ) = 10 * power3( x( 1 ) - x( 2 ) );
    f( 2 ) = power2( tan( x( 2 ) - x( 3 ) ) );
    f( 3 ) = x( 3 ) - 1;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    real_type j0 = 2 * ( exp( x( 0 ) ) - x( 1 ) ) * exp( x( 0 ) );
    real_type j2 = 2 * sin( x( 2 ) - x( 3 ) ) / power3( cos( x( 2 ) - x( 3 ) ) );

    J.insert( 0, 0 ) = j0;
    J.insert( 0, 1 ) = 2 * ( x( 1 ) - exp( x( 0 ) ) );
    J.insert( 0, 2 ) = 0;
    J.insert( 0, 3 ) = 0;

    J.insert( 1, 0 ) = 0;
    J.insert( 1, 1 ) = 30 * power2( x( 1 ) - x( 2 ) );
    J.insert( 1, 2 ) = -30 * power2( x( 1 ) - x( 2 ) );  // Correzione: -30*(x1-x2)^2 invece di -j0
    J.insert( 1, 3 ) = 0;

    J.insert( 2, 0 ) = 0;
    J.insert( 2, 1 ) = 0;
    J.insert( 2, 2 ) = j2;
    J.insert( 2, 3 ) = -j2;

    J.insert( 3, 0 ) = 0;
    J.insert( 3, 1 ) = 0;
    J.insert( 3, 2 ) = 0;
    J.insert( 3, 3 ) = 1;

    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 4, 2, 2, 2;
  }
};
