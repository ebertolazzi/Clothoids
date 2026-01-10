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

class PowellBadlyScaledFunction : public NonlinearSystem
{
  real_type const x0e;
  real_type const x1e;
  real_type const scale;

public:
  PowellBadlyScaledFunction()
    : NonlinearSystem(
        "Powell badly scaled function",
        "@inbook{Powell:1970,\n"
        "  title     = {Numerical methods for nonlinear algebraic "
        "equations},\n"
        "  booktitle = {Proceedings of a {C}onference, {U}niversity of "
        "{E}ssex,\n"
        "              {C}olchester, 6--7 {J}anuary 1969},\n"
        "  chapter   = {An hybrid method for non linear equations},\n"
        "  editor    = {Rabinowitz, Philip},\n"
        "  publisher = {Gordon and Breach Science Publishers, London-New "
        "York-Paris},\n"
        "  year      = {1970},\n"
        "  pages     = {87--114},\n"
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
    //, x0e(1e-6)
    //, x1e(1e6)
    //, scale(1)
    //, x0e(0.109815932969981745568376164563E-4)
    //, x1e(9.10614673986652401094671049032e8)
    //, scale(10000e-18)
    , x0e( 0.109815932969981745568376164563E-4 )
    , x1e( 9.10614673986652401094671049032 )
    , scale( 10000 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = scale * ( ( x( 0 ) * x( 1 ) ) - ( x0e * x1e ) );
    f( 1 ) = ( exp( -x( 1 ) ) - exp( -x1e ) ) + ( exp( -x( 0 ) ) - exp( -x0e ) );
    // f(0) = 10000*(x(0) * x(1)) - 1;
    // f(1) = exp(-x(1)) + exp(-x(0)) - 1.0001;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    J.insert( 0, 0 ) = scale * x( 1 );
    J.insert( 0, 1 ) = scale * x( 0 );

    J.insert( 1, 0 ) = -exp( -x( 0 ) );
    J.insert( 1, 1 ) = -exp( -x( 1 ) );
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << x0e, x1e;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 0, 100;
  }
};
