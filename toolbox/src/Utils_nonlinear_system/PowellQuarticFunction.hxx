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

class PowellQuarticFunction : public NonlinearSystem
{
public:
  PowellQuarticFunction()
    : NonlinearSystem(
        "Powell's quartic function",
        "@article{Colville:1970,\n"
        "  author    = {Colville, A. R.},\n"
        "  title     = {A comparative study of nonlinear programming "
        "codes},\n"
        "  booktitle = {Proceedings of the {P}rinceton {S}ymposium on\n"
        "               {M}athematical {P}rogramming (1967)},\n"
        "  pages     = {487--501},\n"
        "  publisher = {Princeton Univ. Press, Princeton, N.J.},\n"
        "  year      = {1970},\n"
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
        4 )
  {
  }

  virtual void evaluate( Vector const & x_in, Vector & f ) const override
  {
    real_type x = x_in[0];
    real_type y = x_in[1];
    real_type z = x_in[2];
    real_type w = x_in[3];

    // f(0) = 20*w+200*x+4*power3(x-2*y);
    // f(1) = 10*y-10*z-8*power3(x-2*y);
    // f(2) = -10*y+10*z-40*power3(w-z);
    // f(3) = 2*w+20*x+40*power3(w-z);

    real_type t4  = x - 2.0 * y;
    real_type t5  = t4 * t4;
    real_type t6  = t5 * t4;
    real_type t9  = 10.0 * y;
    real_type t10 = 10.0 * z;
    real_type t13 = w - z;
    real_type t14 = t13 * t13;
    real_type t16 = 40.0 * t14 * t13;
    f( 0 )        = 20.0 * w + 200.0 * x + 4.0 * t6;
    f( 1 )        = t9 - t10 - 8.0 * t6;
    f( 2 )        = -t9 + t10 - t16;
    f( 3 )        = 2.0 * w + 20.0 * x + t16;
  }

  virtual void jacobian( Vector const & x_in, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    real_type x = x_in[0];
    real_type y = x_in[1];
    real_type z = x_in[2];
    real_type w = x_in[3];

    real_type t3  = power2( x - 2.0 * y );
    real_type t6  = 24.0 * t3;
    real_type t10 = power2( w - z );
    real_type t11 = 120.0 * t10;

    J.insert( 0, 0 ) = 200.0 + 12.0 * t3;
    J.insert( 0, 1 ) = -t6;
    J.insert( 0, 3 ) = 20.0;

    J.insert( 1, 0 ) = -t6;
    J.insert( 1, 1 ) = 10.0 + 48.0 * t3;
    J.insert( 1, 2 ) = -10.0;

    J.insert( 2, 1 ) = -10.0;
    J.insert( 2, 2 ) = 10.0 + t11;
    J.insert( 2, 3 ) = -t11;

    J.insert( 3, 0 ) = 20.0;
    J.insert( 3, 2 ) = -t11;
    J.insert( 3, 3 ) = 2.0 + t11;

    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.setZero();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 3, -1, 0, 1;
  }
};
