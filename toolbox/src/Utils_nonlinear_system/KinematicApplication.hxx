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

class KinematicApplication : public NonlinearSystem
{
  real_type a1_0, a2_0, a3_0, a4_0, a5_0, a6_0, a7_0, a8_0, a9_0, a10_0, a11_0, a12_0, a13_0, a14_0, a15_0, a16_0,
    a17_0;

  real_type a1_1, a2_1, a3_1, a4_1, a5_1, a6_1, a7_1, a8_1, a9_1, a10_1, a11_1, a12_1, a13_1, a14_1, a15_1, a16_1,
    a17_1;

  real_type a1_2, a2_2, a3_2, a4_2, a5_2, a6_2, a7_2, a8_2, a9_2, a10_2, a11_2, a12_2, a13_2, a14_2, a15_2, a16_2,
    a17_2;

  real_type a1_3, a2_3, a3_3, a4_3, a5_3, a6_3, a7_3, a8_3, a9_3, a10_3, a11_3, a12_3, a13_3, a14_3, a15_3, a16_3,
    a17_3;

public:
  KinematicApplication()
    : NonlinearSystem(
        "Kinematic Application",
        "@article{Morgan:1987,\n"
        "  author  = {Alexander Morgan and Andrew Sommese},\n"
        "  title   = {Computing all solutions to polynomial\n"
        "             systems using homotopy continuation},\n"
        "  journal = {Applied Mathematics and Computation},\n"
        "  volume  = {24},\n"
        "  number  = {2},\n"
        "  pages   = {115--138},\n"
        "  year    = {1987},\n"
        "  doi     = {10.1016/0096-3003(87)90064-6}\n"
        "}\n\n"
        "@article{Hentenryck:1997,\n"
        "  author  = {Van Hentenryck, P. and McAllester, D. "
        "and Kapur, D.},\n"
        "  title   = {Solving Polynomial Systems Using a "
        "Branch and Prune Approach},\n"
        "  journal = {SIAM Journal on Numerical Analysis},\n"
        "  year    = {1997},\n"
        "  volume  = {34},\n"
        "  number  = {2},\n"
        "  pages   = {797-827},\n"
        "  doi = {10.1137/S0036142995281504}\n"
        "}\n",
        8 )
  {
    a1_0  = -0.249150680;
    a2_0  = +1.609135400;
    a3_0  = +0.279423430;
    a4_0  = +1.434801600;
    a5_0  = +0.000000000;
    a6_0  = +0.400263840;
    a7_0  = -0.800527680;
    a8_0  = +0.000000000;
    a9_0  = +0.074052388;
    a10_0 = -0.083050031;
    a11_0 = -0.386159610;
    a12_0 = -0.755266030;
    a13_0 = +0.504201680;
    a14_0 = -1.091628700;
    a15_0 = +0.000000000;
    a16_0 = +0.049207290;
    a17_0 = +0.049207290;

    a1_1  = +0.125016350;
    a2_1  = -0.686607360;
    a3_1  = -0.119228120;
    a4_1  = -0.719940470;
    a5_1  = -0.432419270;
    a6_1  = +0.000000000;
    a7_1  = +0.000000000;
    a8_1  = -0.864838550;
    a9_1  = -0.037157270;
    a10_1 = +0.035436896;
    a11_1 = +0.085383482;
    a12_1 = +0.000000000;
    a13_1 = -0.039251967;
    a14_1 = +0.000000000;
    a15_1 = -0.432419270;
    a16_1 = +0.000000000;
    a17_1 = +0.013873010;

    a1_2  = -0.635550077;
    a2_2  = -0.115719920;
    a3_2  = -0.666404480;
    a4_2  = +0.110362110;
    a5_2  = +0.290702030;
    a6_2  = +1.258776700;
    a7_2  = -0.629388360;
    a8_2  = +0.581404060;
    a9_2  = +0.195946620;
    a10_2 = -1.228034200;
    a11_2 = +0.000000000;
    a12_2 = -0.079034221;
    a13_2 = +0.026387877;
    a14_2 = -0.057131430;
    a15_2 = -1.162808100;
    a16_2 = +1.258776700;
    a17_2 = +2.162575000;

    a1_3  = +1.48947730;
    a2_3  = +0.23062341;
    a3_3  = +1.32810730;
    a4_3  = -0.25864503;
    a5_3  = +1.16517200;
    a6_3  = -0.26908494;
    a7_3  = +0.53816987;
    a8_3  = +0.58258598;
    a9_3  = -0.20816985;
    a10_3 = +2.68683200;
    a11_3 = -0.69910317;
    a12_3 = +0.35744413;
    a13_3 = +1.24991170;
    a14_3 = +1.46773600;
    a15_3 = +1.16517200;
    a16_3 = +1.07633970;
    a17_3 = -0.69686809;
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = x( 0 ) * x( 0 ) + x( 1 ) * x( 1 ) - 1;
    f( 1 ) = x( 1 ) * x( 1 ) + x( 2 ) * x( 2 ) - 1;
    f( 2 ) = x( 2 ) * x( 2 ) + x( 3 ) * x( 3 ) - 1;
    f( 3 ) = x( 3 ) * x( 3 ) + x( 4 ) * x( 4 ) - 1;
    f( 4 ) = a1_0 * x( 0 ) * x( 1 ) + a2_0 * x( 0 ) * x( 3 ) + a3_0 * x( 1 ) * x( 2 ) + a4_0 * x( 1 ) * x( 3 ) +
             a5_0 * x( 1 ) * x( 6 ) + a6_0 * x( 4 ) * x( 7 ) + a7_0 * x( 5 ) * x( 6 ) + a8_0 * x( 5 ) * x( 7 ) +
             a9_0 * x( 0 ) + a10_0 * x( 1 ) + a11_0 * x( 2 ) + a12_0 * x( 3 ) + a13_0 * x( 4 ) + a14_0 * x( 5 ) +
             a15_0 * x( 6 ) + a16_0 * x( 7 ) + a17_0;
    f( 5 ) = a1_1 * x( 0 ) * x( 1 ) + a2_1 * x( 0 ) * x( 3 ) + a3_1 * x( 1 ) * x( 2 ) + a4_1 * x( 1 ) * x( 3 ) +
             a5_1 * x( 1 ) * x( 6 ) + a6_1 * x( 4 ) * x( 7 ) + a7_1 * x( 5 ) * x( 6 ) + a8_1 * x( 5 ) * x( 7 ) +
             a9_1 * x( 0 ) + a10_1 * x( 1 ) + a11_1 * x( 2 ) + a12_1 * x( 3 ) + a13_1 * x( 4 ) + a14_1 * x( 5 ) +
             a15_1 * x( 6 ) + a16_1 * x( 7 ) + a17_1;
    f( 6 ) = a1_2 * x( 0 ) * x( 1 ) + a2_2 * x( 0 ) * x( 3 ) + a3_2 * x( 1 ) * x( 2 ) + a4_2 * x( 1 ) * x( 3 ) +
             a5_2 * x( 1 ) * x( 6 ) + a6_2 * x( 4 ) * x( 7 ) + a7_2 * x( 5 ) * x( 6 ) + a8_2 * x( 5 ) * x( 7 ) +
             a9_2 * x( 0 ) + a10_2 * x( 1 ) + a11_2 * x( 2 ) + a12_2 * x( 3 ) + a13_2 * x( 4 ) + a14_2 * x( 5 ) +
             a15_2 * x( 6 ) + a16_2 * x( 7 ) + a17_2;
    f( 7 ) = a1_3 * x( 0 ) * x( 1 ) + a2_3 * x( 0 ) * x( 3 ) + a3_3 * x( 1 ) * x( 2 ) + a4_3 * x( 1 ) * x( 3 ) +
             a5_3 * x( 1 ) * x( 6 ) + a6_3 * x( 4 ) * x( 7 ) + a7_3 * x( 5 ) * x( 6 ) + a8_3 * x( 5 ) * x( 7 ) +
             a9_3 * x( 0 ) + a10_3 * x( 1 ) + a11_3 * x( 2 ) + a12_3 * x( 3 ) + a13_3 * x( 4 ) + a14_3 * x( 5 ) +
             a15_3 * x( 6 ) + a16_3 * x( 7 ) + a17_3;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    // f(0) = x(0)*x(0) + x(1)*x(1) - 1;
    J.insert( 0, 0 ) = 2 * x( 0 );
    J.insert( 0, 1 ) = 2 * x( 1 );

    // f(1) = x(1)*x(1) + x(2)*x(2) - 1;
    J.insert( 1, 1 ) = 2 * x( 1 );
    J.insert( 1, 2 ) = 2 * x( 2 );

    // f(2) = x(2)*x(2) + x(3)*x(3) - 1;
    J.insert( 2, 2 ) = 2 * x( 2 );
    J.insert( 2, 3 ) = 2 * x( 3 );

    // f(3) = x(3)*x(3) + x(4)*x(4) - 1;
    J.insert( 3, 3 ) = 2 * x( 3 );
    J.insert( 3, 4 ) = 2 * x( 4 );

    J.insert( 4, 0 ) = a1_0 * x( 1 ) + a2_0 * x( 3 ) + a9_0;
    J.insert( 4, 1 ) = a1_0 * x( 0 ) + a3_0 * x( 2 ) + a4_0 * x( 3 ) + a5_0 * x( 6 ) + a10_0;
    J.insert( 4, 2 ) = a3_0 * x( 1 ) + a11_0;
    J.insert( 4, 3 ) = a2_0 * x( 0 ) + a4_0 * x( 1 ) + a12_0;
    J.insert( 4, 4 ) = a6_0 * x( 7 ) + a13_0;
    J.insert( 4, 5 ) = a7_0 * x( 6 ) + a8_0 * x( 7 ) + a14_0;
    J.insert( 4, 6 ) = a5_0 * x( 1 ) + a7_0 * x( 5 ) + a15_0;
    J.insert( 4, 7 ) = a6_0 * x( 4 ) + a8_0 * x( 5 ) + a16_0;

    J.insert( 5, 0 ) = a1_1 * x( 1 ) + a2_1 * x( 3 ) + a9_1;
    J.insert( 5, 1 ) = a1_1 * x( 0 ) + a3_1 * x( 2 ) + a4_1 * x( 3 ) + a5_1 * x( 6 ) + a10_1;
    J.insert( 5, 2 ) = a3_1 * x( 1 ) + a11_1;
    J.insert( 5, 3 ) = a2_1 * x( 0 ) + a4_1 * x( 1 ) + a12_1;
    J.insert( 5, 4 ) = a6_1 * x( 7 ) + a13_1;
    J.insert( 5, 5 ) = a7_1 * x( 6 ) + a8_1 * x( 7 ) + a14_1;
    J.insert( 5, 6 ) = a5_1 * x( 1 ) + a7_1 * x( 6 ) + a15_1;
    J.insert( 5, 7 ) = a6_1 * x( 4 ) + a8_1 * x( 5 ) + a16_1;

    J.insert( 6, 0 ) = a1_2 * x( 1 ) + a2_2 * x( 3 ) + a9_2;
    J.insert( 6, 1 ) = a1_2 * x( 0 ) + a3_2 * x( 2 ) + a4_2 * x( 3 ) + a5_2 * x( 6 ) + a10_2;
    J.insert( 6, 2 ) = a3_2 * x( 1 ) + a11_2;
    J.insert( 6, 3 ) = a2_2 * x( 0 ) + a4_2 * x( 1 ) + a12_2;
    J.insert( 6, 4 ) = a6_2 * x( 7 ) + a13_2;
    J.insert( 6, 5 ) = a7_2 * x( 6 ) + a8_2 * x( 7 ) + a14_2;
    J.insert( 6, 6 ) = a5_2 * x( 1 ) + a7_2 * x( 5 ) + a15_2;
    J.insert( 6, 7 ) = a6_2 * x( 4 ) + a8_2 * x( 5 ) + a16_2;

    J.insert( 7, 0 ) = a1_3 * x( 1 ) + a2_3 * x( 3 ) + a9_3;
    J.insert( 7, 1 ) = a1_3 * x( 0 ) + a3_3 * x( 2 ) + a4_3 * x( 3 ) + a5_3 * x( 6 ) + a10_3;
    J.insert( 7, 2 ) = a3_3 * x( 1 ) + a11_3;
    J.insert( 7, 3 ) = a2_3 * x( 0 ) + a4_3 * x( 1 ) + a12_3;
    J.insert( 7, 4 ) = a6_3 * x( 7 ) + a13_3;
    J.insert( 7, 5 ) = a7_3 * x( 6 ) + a8_3 * x( 7 ) + a14_3;
    J.insert( 7, 6 ) = a5_3 * x( 1 ) + a7_3 * x( 5 ) + a15_3;
    J.insert( 7, 7 ) = a6_3 * x( 4 ) + a8_3 * x( 5 ) + a16_3;

    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 0.5 );
  }

  virtual void check_if_admissible( Vector const & x ) const override
  {
    for ( integer i = 0; i < n; ++i ) UTILS_ASSERT( std::abs( x( i ) ) < 1000, "Bad range" );
  }

  virtual void bounding_box( Vector & L, Vector & U ) const override
  {
    U.fill( 1000 );
    L.fill( -10000 );
  }
};
