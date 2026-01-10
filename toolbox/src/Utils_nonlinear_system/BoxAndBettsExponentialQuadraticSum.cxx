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

class BoxAndBettsExponentialQuadraticSum : public NonlinearSystem
{
public:
  BoxAndBettsExponentialQuadraticSum()
    : NonlinearSystem(
        "Box and Betts exponential quadratic sum",
        "@article{Box:1966,\n"
        "  author  = {Box, M. J.},\n"
        "  title   = {A Comparison of Several Current Optimization "
        "Methods,\n"
        "             and the use of Transformations in Constrained "
        "Problems},\n"
        "  journal = {The Computer Journal},\n"
        "  volume  = {9},\n"
        "  number  = {1},\n"
        "  pages   = {67--77},\n"
        "  year    = {1966},\n"
        "  doi     = {10.1093/comjnl/9.1.67},\n"
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
        3 )
  {
  }

  virtual void evaluate( Vector const & x_in, Vector & f ) const override
  {
    real_type const & x    = x_in( 0 );
    real_type const & y    = x_in( 1 );
    real_type const & z    = x_in( 2 );
    real_type         t2   = exp( -x / 10.0 );
    real_type         t4   = exp( -y / 10.0 );
    real_type         t5   = exp( -1.0 / 10.0 );
    real_type         t6   = exp( -1.0 );
    real_type         t7   = t5 - t6;
    real_type         t9   = t2 - t4 - t7 * z;
    real_type         t13  = exp( -x / 5.0 );
    real_type         t15  = exp( -y / 5.0 );
    real_type         t16  = exp( -1.0 / 5.0 );
    real_type         t17  = exp( -2.0 );
    real_type         t18  = t16 - t17;
    real_type         t20  = t13 - t15 - t18 * z;
    real_type         t24  = exp( -3.0 / 10.0 * x );
    real_type         t26  = exp( -3.0 / 10.0 * y );
    real_type         t27  = exp( -3.0 / 10.0 );
    real_type         t28  = exp( -3.0 );
    real_type         t29  = t27 - t28;
    real_type         t31  = t24 - t26 - t29 * z;
    real_type         t35  = exp( -2.0 / 5.0 * x );
    real_type         t37  = exp( -2.0 / 5.0 * y );
    real_type         t38  = exp( -2.0 / 5.0 );
    real_type         t39  = exp( -4.0 );
    real_type         t40  = t38 - t39;
    real_type         t42  = t35 - t37 - t40 * z;
    real_type         t46  = exp( -x / 2.0 );
    real_type         t48  = exp( -y / 2.0 );
    real_type         t49  = exp( -1.0 / 2.0 );
    real_type         t50  = exp( -5.0 );
    real_type         t51  = t49 - t50;
    real_type         t53  = t46 - t48 - t51 * z;
    real_type         t56  = exp( -3.0 / 5.0 * x );
    real_type         t58  = exp( -3.0 / 5.0 * y );
    real_type         t59  = exp( -3.0 / 5.0 );
    real_type         t60  = exp( -6.0 );
    real_type         t61  = t59 - t60;
    real_type         t63  = t56 - t58 - t61 * z;
    real_type         t67  = exp( -7.0 / 10.0 * x );
    real_type         t69  = exp( -7.0 / 10.0 * y );
    real_type         t70  = exp( -7.0 / 10.0 );
    real_type         t71  = exp( -7.0 );
    real_type         t72  = t70 - t71;
    real_type         t74  = t67 - t69 - t72 * z;
    real_type         t78  = exp( -4.0 / 5.0 * x );
    real_type         t80  = exp( -4.0 / 5.0 * y );
    real_type         t81  = exp( -4.0 / 5.0 );
    real_type         t82  = exp( -8.0 );
    real_type         t83  = t81 - t82;
    real_type         t85  = t78 - t80 - t83 * z;
    real_type         t89  = exp( -9.0 / 10.0 * x );
    real_type         t91  = exp( -9.0 / 10.0 * y );
    real_type         t92  = exp( -9.0 / 10.0 );
    real_type         t93  = exp( -9.0 );
    real_type         t94  = t92 - t93;
    real_type         t96  = t89 - t91 - t94 * z;
    real_type         t99  = exp( -x );
    real_type         t100 = exp( -y );
    real_type         t101 = exp( -10.0 );
    real_type         t102 = t6 - t101;
    real_type         t104 = t99 - t100 - t102 * z;
    f( 0 ) = -t9 * t2 / 5.0 - 2.0 / 5.0 * t20 * t13 - 3.0 / 5.0 * t31 * t24 - 4.0 / 5.0 * t42 * t35 - t53 * t46 -
             6.0 / 5.0 * t63 * t56 - 7.0 / 5.0 * t74 * t67 - 8.0 / 5.0 * t85 * t78 - 9.0 / 5.0 * t96 * t89 -
             2.0 * t104 * t99;
    f( 1 ) = t9 * t4 / 5.0 + 2.0 / 5.0 * t20 * t15 + 3.0 / 5.0 * t31 * t26 + 4.0 / 5.0 * t42 * t37 + t53 * t48 +
             6.0 / 5.0 * t63 * t58 + 7.0 / 5.0 * t74 * t69 + 8.0 / 5.0 * t85 * t80 + 9.0 / 5.0 * t96 * t91 +
             2.0 * t104 * t100;
    f( 2 ) = -2.0 * t9 * t7 - 2.0 * t20 * t18 - 2.0 * t31 * t29 - 2.0 * t42 * t40 - 2.0 * t53 * t51 - 2.0 * t63 * t61 -
             2.0 * t74 * t72 - 2.0 * t85 * t83 - 2.0 * t96 * t94 - 2.0 * t104 * t102;
  }

  virtual void jacobian( Vector const & x_in, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    real_type const & x    = x_in( 0 );
    real_type const & y    = x_in( 1 );
    real_type const & z    = x_in( 2 );
    real_type         t2   = exp( -x / 10.0 );
    real_type         t3   = t2 * t2;
    real_type         t6   = exp( -y / 10.0 );
    real_type         t7   = exp( -1.0 / 10.0 );
    real_type         t8   = exp( -1.0 );
    real_type         t9   = t7 - t8;
    real_type         t11  = t2 - t6 - t9 * z;
    real_type         t15  = exp( -x / 5.0 );
    real_type         t16  = t15 * t15;
    real_type         t19  = exp( -y / 5.0 );
    real_type         t20  = exp( -1.0 / 5.0 );
    real_type         t21  = exp( -2.0 );
    real_type         t22  = t20 - t21;
    real_type         t24  = t15 - t19 - t22 * z;
    real_type         t28  = exp( -3.0 / 10.0 * x );
    real_type         t29  = t28 * t28;
    real_type         t32  = exp( -3.0 / 10.0 * y );
    real_type         t33  = exp( -3.0 / 10.0 );
    real_type         t34  = exp( -3.0 );
    real_type         t35  = t33 - t34;
    real_type         t37  = t28 - t32 - t35 * z;
    real_type         t41  = exp( -2.0 / 5.0 * x );
    real_type         t42  = t41 * t41;
    real_type         t45  = exp( -2.0 / 5.0 * y );
    real_type         t46  = exp( -2.0 / 5.0 );
    real_type         t47  = exp( -4.0 );
    real_type         t48  = t46 - t47;
    real_type         t50  = t41 - t45 - t48 * z;
    real_type         t54  = exp( -x / 2.0 );
    real_type         t55  = t54 * t54;
    real_type         t58  = exp( -y / 2.0 );
    real_type         t59  = exp( -1.0 / 2.0 );
    real_type         t60  = exp( -5.0 );
    real_type         t61  = t59 - t60;
    real_type         t63  = t54 - t58 - t61 * z;
    real_type         t67  = exp( -3.0 / 5.0 * x );
    real_type         t68  = t67 * t67;
    real_type         t71  = exp( -3.0 / 5.0 * y );
    real_type         t72  = exp( -3.0 / 5.0 );
    real_type         t73  = exp( -6.0 );
    real_type         t74  = t72 - t73;
    real_type         t76  = t67 - t71 - t74 * z;
    real_type         t80  = exp( -7.0 / 10.0 * x );
    real_type         t81  = t80 * t80;
    real_type         t84  = exp( -7.0 / 10.0 * y );
    real_type         t85  = exp( -7.0 / 10.0 );
    real_type         t86  = exp( -7.0 );
    real_type         t87  = t85 - t86;
    real_type         t89  = t80 - t84 - t87 * z;
    real_type         t93  = exp( -4.0 / 5.0 * x );
    real_type         t94  = t93 * t93;
    real_type         t97  = exp( -4.0 / 5.0 * y );
    real_type         t98  = exp( -4.0 / 5.0 );
    real_type         t99  = exp( -8.0 );
    real_type         t100 = t98 - t99;
    real_type         t102 = t93 - t97 - t100 * z;
    real_type         t106 = exp( -9.0 / 10.0 * x );
    real_type         t107 = t106 * t106;
    real_type         t110 = exp( -9.0 / 10.0 * y );
    real_type         t111 = exp( -9.0 / 10.0 );
    real_type         t112 = exp( -9.0 );
    real_type         t113 = t111 - t112;
    real_type         t115 = t106 - t110 - t113 * z;
    real_type         t118 = exp( -x );
    real_type         t119 = t118 * t118;
    real_type         t121 = exp( -y );
    real_type         t122 = exp( -10.0 );
    real_type         t123 = t8 - t122;
    real_type         t125 = t118 - t121 - t123 * z;
    real_type t128 = t3 / 50.0 + t11 * t2 / 50.0 + 2.0 / 25.0 * t16 + 2.0 / 25.0 * t24 * t15 + 9.0 / 50.0 * t29 +
                     9.0 / 50.0 * t37 * t28 + 8.0 / 25.0 * t42 + 8.0 / 25.0 * t50 * t41 + t55 / 2.0 + t63 * t54 / 2.0 +
                     18.0 / 25.0 * t68 + 18.0 / 25.0 * t76 * t67 + 49.0 / 50.0 * t81 + 49.0 / 50.0 * t89 * t80 +
                     32.0 / 25.0 * t94 + 32.0 / 25.0 * t102 * t93 + 81.0 / 50.0 * t107 + 81.0 / 50.0 * t115 * t106 +
                     2.0 * t119 + 2.0 * t125 * t118;
    real_type t149 = -t6 * t2 / 50.0 - 2.0 / 25.0 * t19 * t15 - 9.0 / 50.0 * t32 * t28 - 8.0 / 25.0 * t45 * t41 -
                     t58 * t54 / 2.0 - 18.0 / 25.0 * t71 * t67 - 49.0 / 50.0 * t84 * t80 - 32.0 / 25.0 * t97 * t93 -
                     81.0 / 50.0 * t110 * t106 - 2.0 * t121 * t118;
    real_type t169 = t9 * t2 / 5.0 + 2.0 / 5.0 * t15 * t22 + 3.0 / 5.0 * t28 * t35 + 4.0 / 5.0 * t41 * t48 + t54 * t61 +
                     6.0 / 5.0 * t74 * t67 + 7.0 / 5.0 * t80 * t87 + 8.0 / 5.0 * t93 * t100 + 9.0 / 5.0 * t106 * t113 +
                     2.0 * t118 * t123;
    real_type t170 = t6 * t6;
    real_type t174 = t19 * t19;
    real_type t178 = t32 * t32;
    real_type t182 = t45 * t45;
    real_type t186 = t58 * t58;
    real_type t190 = t71 * t71;
    real_type t194 = t84 * t84;
    real_type t198 = t97 * t97;
    real_type t202 = t110 * t110;
    real_type t206 = t121 * t121;
    real_type t210 = t170 / 50.0 - t11 * t6 / 50.0 + 2.0 / 25.0 * t174 - 2.0 / 25.0 * t24 * t19 + 9.0 / 50.0 * t178 -
                     9.0 / 50.0 * t37 * t32 + 8.0 / 25.0 * t182 - 8.0 / 25.0 * t50 * t45 + t186 / 2.0 -
                     t63 * t58 / 2.0 + 18.0 / 25.0 * t190 - 18.0 / 25.0 * t76 * t71 + 49.0 / 50.0 * t194 -
                     49.0 / 50.0 * t89 * t84 + 32.0 / 25.0 * t198 - 32.0 / 25.0 * t102 * t97 + 81.0 / 50.0 * t202 -
                     81.0 / 50.0 * t115 * t110 + 2.0 * t206 - 2.0 * t125 * t121;
    real_type t230 = -t6 * t9 / 5.0 - 2.0 / 5.0 * t19 * t22 - 3.0 / 5.0 * t32 * t35 - 4.0 / 5.0 * t45 * t48 -
                     t58 * t61 - 6.0 / 5.0 * t71 * t74 - 7.0 / 5.0 * t84 * t87 - 8.0 / 5.0 * t97 * t100 -
                     9.0 / 5.0 * t110 * t113 - 2.0 * t121 * t123;
    real_type t231   = t9 * t9;
    real_type t232   = t22 * t22;
    real_type t233   = t35 * t35;
    real_type t234   = t48 * t48;
    real_type t235   = t61 * t61;
    real_type t236   = t74 * t74;
    real_type t237   = t87 * t87;
    real_type t238   = t100 * t100;
    real_type t239   = t113 * t113;
    real_type t240   = t123 * t123;
    J.insert( 0, 0 ) = t128;
    J.insert( 0, 1 ) = t149;
    J.insert( 0, 2 ) = t169;
    J.insert( 1, 0 ) = t149;
    J.insert( 1, 1 ) = t210;
    J.insert( 1, 2 ) = t230;
    J.insert( 2, 0 ) = t169;
    J.insert( 2, 1 ) = t230;
    J.insert( 2, 2 ) = 2.0 * t231 + 2.0 * t232 + 2.0 * t233 + 2.0 * t234 + 2.0 * t235 + 2.0 * t236 + 2.0 * t237 +
                       2.0 * t238 + 2.0 * t239 + 2.0 * t240;
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
    x0 << 0, 10, 20;
  }
};
