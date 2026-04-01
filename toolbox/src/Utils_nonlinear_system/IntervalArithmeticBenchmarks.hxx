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

class IntervalArithmeticBenchmarks : public NonlinearSystem
{
  real_type a1, a2, a3, a4, a5, a6, a7, a8, a9, a10;
  real_type b1, b2, b3, b4, b5, b6, b7, b8, b9, b10;

public:
  IntervalArithmeticBenchmarks()
    : NonlinearSystem(
        "Interval Arithmetic Benchmarks",
        "@article{Morgan:1987,\n"
        "  author  = {Alexander Morgan and Andrew Sommese},\n"
        "  title   = {Computing all solutions to polynomial\n"
        "             systems using homotopy continuation},\n"
        "  journal = {Applied Mathematics and Computation},\n"
        "  volume  = {24},\n"
        "  number  = {2},\n"
        "  pages   = {115--138},\n"
        "  year    = {1987},\n"
        "  issn    = {0096-3003},\n"
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
        10 )
  {
    a1  = 0.25428722;
    a2  = 0.37842197;
    a3  = 0.27162577;
    a4  = 0.19807914;
    a5  = 0.44166728;
    a6  = 0.14654113;
    a7  = 0.42937161;
    a8  = 0.07056438;
    a9  = 0.34504906;
    a10 = 0.42651102;

    b1  = 0.18324757;
    b2  = 0.16275449;
    b3  = 0.16955071;
    b4  = 0.15585316;
    b5  = 0.19950920;
    b6  = 0.18922793;
    b7  = 0.21180486;
    b8  = 0.17081208;
    b9  = 0.19612740;
    b10 = 0.21466544;
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = x( 0 ) - a1 - b1 * x( 3 ) * x( 2 ) * x( 8 );
    f( 1 ) = x( 1 ) - a2 - b2 * x( 0 ) * x( 9 ) * x( 5 );
    f( 2 ) = x( 2 ) - a3 - b3 * x( 0 ) * x( 1 ) * x( 9 );
    f( 3 ) = x( 3 ) - a4 - b4 * x( 6 ) * x( 0 ) * x( 5 );
    f( 4 ) = x( 4 ) - a5 - b5 * x( 6 ) * x( 5 ) * x( 2 );
    f( 5 ) = x( 5 ) - a6 - b6 * x( 7 ) * x( 4 ) * x( 9 );
    f( 6 ) = x( 6 ) - a7 - b7 * x( 1 ) * x( 4 ) * x( 7 );
    f( 7 ) = x( 7 ) - a8 - b8 * x( 0 ) * x( 6 ) * x( 5 );
    f( 8 ) = x( 8 ) - a9 - b9 * x( 9 ) * x( 5 ) * x( 7 );
    f( 9 ) = x( 9 ) - a10 - b10 * x( 3 ) * x( 7 ) * x( 0 );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    J.insert( 0, 0 ) = 1;
    J.insert( 0, 2 ) = -b1 * x( 3 ) * x( 8 );
    J.insert( 0, 3 ) = -b1 * x( 2 ) * x( 8 );
    J.insert( 0, 8 ) = -b1 * x( 3 ) * x( 2 );

    J.insert( 1, 1 ) = 1;
    J.insert( 1, 0 ) = -b2 * x( 9 ) * x( 5 );
    J.insert( 1, 5 ) = -b2 * x( 0 ) * x( 9 );
    J.insert( 1, 9 ) = -b2 * x( 0 ) * x( 5 );

    J.insert( 2, 2 ) = 1;
    J.insert( 2, 0 ) = -b3 * x( 1 ) * x( 9 );
    J.insert( 2, 1 ) = -b3 * x( 0 ) * x( 9 );
    J.insert( 2, 9 ) = -b3 * x( 0 ) * x( 1 );

    J.insert( 3, 3 ) = 1;
    J.insert( 3, 0 ) = -b4 * x( 6 ) * x( 5 );
    J.insert( 3, 5 ) = -b4 * x( 6 ) * x( 0 );
    J.insert( 3, 6 ) = -b4 * x( 0 ) * x( 5 );

    J.insert( 4, 4 ) = 1;
    J.insert( 4, 2 ) = -b5 * x( 6 ) * x( 5 );
    J.insert( 4, 5 ) = -b5 * x( 6 ) * x( 2 );
    J.insert( 4, 6 ) = -b5 * x( 5 ) * x( 2 );

    J.insert( 5, 5 ) = 1;
    J.insert( 5, 4 ) = -b6 * x( 7 ) * x( 9 );
    J.insert( 5, 7 ) = -b6 * x( 4 ) * x( 9 );
    J.insert( 5, 9 ) = -b6 * x( 7 ) * x( 4 );

    J.insert( 6, 6 ) = 1;
    J.insert( 6, 1 ) = -b7 * x( 4 ) * x( 7 );
    J.insert( 6, 4 ) = -b7 * x( 1 ) * x( 7 );
    J.insert( 6, 7 ) = -b7 * x( 1 ) * x( 4 );

    J.insert( 7, 7 ) = 1;
    J.insert( 7, 0 ) = -b8 * x( 6 ) * x( 5 );
    J.insert( 7, 6 ) = -b8 * x( 0 ) * x( 5 );
    J.insert( 7, 5 ) = -b8 * x( 0 ) * x( 6 );

    J.insert( 8, 8 ) = 1;
    J.insert( 8, 5 ) = -b9 * x( 9 ) * x( 7 );
    J.insert( 8, 7 ) = -b9 * x( 9 ) * x( 5 );
    J.insert( 8, 9 ) = -b9 * x( 5 ) * x( 7 );

    J.insert( 9, 9 ) = 1;
    J.insert( 9, 3 ) = -b10 * x( 7 ) * x( 0 );
    J.insert( 9, 7 ) = -b10 * x( 3 ) * x( 0 );
    J.insert( 9, 0 ) = -b10 * x( 3 ) * x( 7 );

    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.setZero();
  }
};
