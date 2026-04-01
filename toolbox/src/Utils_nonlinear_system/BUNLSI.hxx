
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

#define BUNLSI_BIBTEX                                                 \
  "@article{Buzzi:1986,\n"                                            \
  "  author  = {Guido Buzzi Ferraris and Enrico Tronconi},\n"         \
  "  title   = {Bunlsi—{A} fortran program for solution of systems\n" \
  "             of nonlinear algebraic equations},\n"                 \
  "  journal = {Computers \\& Chemical Engineering},\n"               \
  "  volume  = {10},\n"                                               \
  "  number  = {2},\n"                                                \
  "  pages   = {129--141},\n"                                         \
  "  year    = {1986},\n"                                             \
  "  doi     = {10.1016/0098-1354(86)85025-6},\n"                     \
  "}\n"

class BUNLSI5 : public NonlinearSystem
{
public:
  BUNLSI5() : NonlinearSystem( "BUNLSI example 5", BUNLSI_BIBTEX, 8 ) {}

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type x1 = x( 0 );
    real_type x2 = x( 1 );
    real_type x3 = x( 2 );
    real_type x4 = x( 3 );
    real_type x5 = x( 4 );
    real_type x6 = x( 5 );
    real_type x7 = x( 6 );
    real_type x8 = x( 7 );
    f( 0 )       = x1 - 1;
    f( 1 )       = x2 - sqrt( x1 ) - exp( x1 ) - 15;
    f( 2 )       = x3 - x1 * x2 / 100 - sin( x1 ) - 1;
    f( 3 )       = x4 - power2( x2 + x1 + x3 / 2 ) + 150;
    f( 4 )       = x5 - ( x4 - x2 ) / ( x1 * x2 * x3 );
    f( 5 )       = x6 - x5 * pow( x1, 1.0 / 3.0 ) - exp( x5 );
    f( 6 )       = x7 - ( x1 - sqrt( x4 ) - x5 * x5 ) * exp( x5 );
    f( 7 )       = x8 - x1 - x5 - x6 - x3;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    real_type x1 = x( 0 );
    real_type x2 = x( 1 );
    real_type x3 = x( 2 );
    real_type x4 = x( 3 );
    real_type x5 = x( 4 );
    // real_type x6 = x(5);
    // real_type x7 = x(6);
    // real_type x8 = x(7);

    J.insert( 0, 0 ) = 1;

    J.insert( 1, 0 ) = -0.5 / sqrt( x1 ) - exp( x1 );
    J.insert( 1, 1 ) = 1;

    J.insert( 2, 0 ) = -x2 / 100 - cos( x1 );
    J.insert( 2, 1 ) = -x1 / 100;
    J.insert( 2, 2 ) = 1;

    J.insert( 3, 0 ) = -2 * ( x2 + x1 + x3 / 2 );
    J.insert( 3, 1 ) = -2 * ( x2 + x1 + x3 / 2 );
    J.insert( 3, 2 ) = -( x2 + x1 + x3 / 2 );
    J.insert( 3, 3 ) = 1;

    J.insert( 4, 0 ) = ( x4 - x2 ) / ( x1 * x1 * x2 * x3 );
    J.insert( 4, 1 ) = x4 / ( x1 * x2 * x2 * x3 );
    J.insert( 4, 2 ) = ( x4 - x2 ) / ( x1 * x2 * x3 * x3 );
    J.insert( 4, 3 ) = -1 / ( x1 * x2 * x3 );
    J.insert( 4, 4 ) = 1;

    J.insert( 5, 0 ) = -x5 * pow( x1, -2.0 / 3.0 ) / 3;
    J.insert( 5, 4 ) = -pow( x1, 1.0 / 3.0 ) - exp( x5 );
    J.insert( 5, 5 ) = 1;

    J.insert( 6, 0 ) = -exp( x5 );
    J.insert( 6, 3 ) = 0.5 * exp( x5 ) / sqrt( x4 );
    J.insert( 6, 4 ) = ( x5 * ( x5 + 2 ) + sqrt( x4 ) - x1 ) * exp( x5 );
    J.insert( 6, 6 ) = 1;

    J.insert( 7, 0 ) = -1;
    J.insert( 7, 2 ) = -1;
    J.insert( 7, 4 ) = -1;
    J.insert( 7, 5 ) = -1;
    J.insert( 7, 7 ) = 1;

    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 1.000000000000000000000000000000000000000, 18.71828182845904523536028747135266249776,
      2.028653803092486959006105196343825624600, 279.8410647514915106118923581352729397598,
      6.876553786347144986146293257311354434906, 976.1568043532779061610180406939220907971,
      -61079.62412297820687494549945165294788674, 986.0620119427175381061704391475772708566;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 1, 20, 2.2, 100, 2, 8, -60, 15;
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class BUNLSI6 : public NonlinearSystem
{
public:
  BUNLSI6() : NonlinearSystem( "BUNLSI example 6", BUNLSI_BIBTEX, 30 ) {}

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type x1  = x( 0 );
    real_type x2  = x( 1 );
    real_type x3  = x( 2 );
    real_type x4  = x( 3 );
    real_type x5  = x( 4 );
    real_type x6  = x( 5 );
    real_type x7  = x( 6 );
    real_type x8  = x( 7 );
    real_type x9  = x( 8 );
    real_type x10 = x( 9 );
    real_type x11 = x( 10 );
    real_type x12 = x( 11 );
    real_type x13 = x( 12 );
    real_type x14 = x( 13 );
    real_type x15 = x( 14 );
    real_type x16 = x( 15 );
    real_type x17 = x( 16 );
    real_type x18 = x( 17 );
    real_type x19 = x( 18 );
    real_type x20 = x( 19 );
    real_type x21 = x( 20 );
    real_type x22 = x( 21 );
    real_type x23 = x( 22 );
    real_type x24 = x( 23 );
    real_type x25 = x( 24 );
    real_type x26 = x( 25 );
    real_type x27 = x( 26 );
    real_type x28 = x( 27 );
    real_type x29 = x( 28 );
    real_type x30 = x( 29 );
    f( 0 )        = x1 - 1;
    f( 1 )        = 6 * x1 - x2;
    f( 2 )        = 5.4 * x1 + x2 - x3;
    f( 3 )        = x2 + x3 - x4;
    f( 4 )        = 100 * ( x1 + x2 - x5 );
    f( 5 )        = x5 + 100 * x4 - x6;
    f( 6 )        = x1 - 0.1 * x4 - x7;
    f( 7 )        = x5 + 90 * x2 - x8;
    f( 8 )        = x1 + x7 - x9;
    f( 9 )        = x9 + x3 - x10;
    f( 10 )       = x1 * x2 - x11;
    f( 11 )       = x11 / x3 - x12;
    f( 12 )       = sqrt( x5 ) + x4 - x13;
    f( 13 )       = log( x1 * x6 ) + x10 - x14;
    f( 14 )       = -sin( x2 ) + log( x2 + x10 ) - x15;
    f( 15 )       = x3 * x9 * x10 - x14 - x16;
    f( 16 )       = 69.1 * x1 - 0.01 * x2 * x5 - x17;
    f( 17 )       = x6 / ( x5 * x15 ) - x18;
    f( 18 )       = sqrt( x10 * x16 ) - x19;
    f( 19 )       = pow( x1, 2.6 ) + x12 * x12 - x20;
    f( 20 )       = x11 * exp( x20 ) - x21;
    f( 21 )       = x5 / x10 + x13 - x22;
    f( 22 )       = 0.1 * x2 + x1 * x10 / x4 - x23;
    f( 23 )       = log( x5 * x8 ) + x2 - x24;
    f( 24 )       = x11 * x24 - x22 - x25;
    f( 25 )       = power2( x11 * x12 * x14 ) - x26 - 1500;
    f( 26 )       = x1 * x5 / 10000 - x27;
    f( 27 )       = power2( x6 - x8 ) / x22 - x28;
    f( 28 )       = x16 * x19 / x21 - x29;
    f( 29 )       = x27 * x28 + x22 + sqrt( x28 ) - x30;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type x1 = x( 0 );
    real_type x2 = x( 1 );
    real_type x3 = x( 2 );
    real_type x4 = x( 3 );
    real_type x5 = x( 4 );
    real_type x6 = x( 5 );
    // real_type x7 = x(6);
    real_type x8  = x( 7 );
    real_type x9  = x( 8 );
    real_type x10 = x( 9 );
    real_type x11 = x( 10 );
    real_type x12 = x( 11 );
    // real_type x13 = x(12);
    real_type x14 = x( 13 );
    real_type x15 = x( 14 );
    real_type x16 = x( 15 );
    // real_type x17 = x(16);
    // real_type x18 = x(17);
    real_type x19 = x( 18 );
    real_type x20 = x( 19 );
    real_type x21 = x( 20 );
    real_type x22 = x( 21 );
    // real_type x23 = x(22);
    real_type x24 = x( 23 );
    // real_type x25 = x(24);
    // real_type x26 = x(25);
    real_type x27 = x( 26 );
    real_type x28 = x( 27 );
    // real_type x29 = x(28);
    // real_type x30 = x(29);

    J.insert( 0, 0 ) = 1;

    J.insert( 1, 0 ) = 6;
    J.insert( 1, 1 ) = -1;

    J.insert( 2, 0 ) = 5.4;
    J.insert( 2, 1 ) = 1;
    J.insert( 2, 2 ) = -1;

    J.insert( 3, 1 ) = 1;
    J.insert( 3, 2 ) = 1;
    J.insert( 3, 3 ) = -1;

    J.insert( 4, 0 ) = 100;
    J.insert( 4, 1 ) = 100;
    J.insert( 4, 4 ) = -100;

    J.insert( 5, 3 ) = 100;
    J.insert( 5, 4 ) = 1;
    J.insert( 5, 5 ) = -1;

    J.insert( 6, 0 ) = 1;
    J.insert( 6, 3 ) = -0.1;
    J.insert( 6, 6 ) = -1;

    J.insert( 7, 1 ) = 90;
    J.insert( 7, 4 ) = 1;
    J.insert( 7, 7 ) = -1;

    J.insert( 8, 0 ) = 1;
    J.insert( 8, 6 ) = 1;
    J.insert( 8, 8 ) = -1;

    J.insert( 9, 2 ) = 1;
    J.insert( 9, 8 ) = 1;
    J.insert( 9, 9 ) = -1;

    J.insert( 10, 0 )  = x2;
    J.insert( 10, 1 )  = x1;
    J.insert( 10, 10 ) = -1;

    J.insert( 11, 2 )  = -x11 / ( x3 * x3 );
    J.insert( 11, 10 ) = 1 / x3;
    J.insert( 11, 11 ) = -1;

    J.insert( 12, 3 )  = 1;
    J.insert( 12, 4 )  = 0.5 / sqrt( x5 );
    J.insert( 12, 12 ) = -1;

    J.insert( 13, 0 )  = 1 / x1;
    J.insert( 13, 5 )  = 1 / x6;
    J.insert( 13, 9 )  = 1;
    J.insert( 13, 13 ) = -1;

    J.insert( 14, 1 )  = 1 / ( x2 + x10 ) - cos( x2 );
    J.insert( 14, 9 )  = 1 / ( x2 + x10 );
    J.insert( 14, 14 ) = -1;

    J.insert( 15, 2 )  = x9 * x10;
    J.insert( 15, 8 )  = x3 * x10;
    J.insert( 15, 9 )  = x3 * x9;
    J.insert( 15, 13 ) = -1;
    J.insert( 15, 15 ) = -1;

    J.insert( 16, 0 )  = 69.1;
    J.insert( 16, 1 )  = -0.01 * x5;
    J.insert( 16, 4 )  = -0.01 * x2;
    J.insert( 16, 16 ) = -1;

    J.insert( 17, 4 )  = -x6 / ( x5 * x5 * x15 );
    J.insert( 17, 5 )  = 1 / ( x5 * x15 );
    J.insert( 17, 14 ) = -x6 / ( x5 * x15 * x15 );
    J.insert( 17, 17 ) = -1;

    // CORREZIONE: indice 17 -> 15 per la derivata rispetto a x16
    J.insert( 18, 9 )  = 0.5 * x16 / sqrt( x10 * x16 );
    J.insert( 18, 15 ) = 0.5 * x10 / sqrt( x10 * x16 );  // Era 17, corretto a 15
    J.insert( 18, 18 ) = -1;

    J.insert( 19, 0 )  = 2.6 * pow( x1, 1.6 );
    J.insert( 19, 11 ) = 2 * x12;
    J.insert( 19, 19 ) = -1;

    J.insert( 20, 10 ) = exp( x20 );
    J.insert( 20, 19 ) = x11 * exp( x20 );
    J.insert( 20, 20 ) = -1;

    J.insert( 21, 4 )  = 1 / x10;
    J.insert( 21, 9 )  = -x5 / ( x10 * x10 );
    J.insert( 21, 12 ) = 1;
    J.insert( 21, 21 ) = -1;

    J.insert( 22, 0 )  = x10 / x4;
    J.insert( 22, 1 )  = 0.1;
    J.insert( 22, 3 )  = -x1 * x10 / ( x4 * x4 );
    J.insert( 22, 9 )  = x1 / x4;
    J.insert( 22, 22 ) = -1;

    J.insert( 23, 1 )  = 1;
    J.insert( 23, 4 )  = 1 / x5;
    J.insert( 23, 7 )  = 1 / x8;
    J.insert( 23, 23 ) = -1;

    J.insert( 24, 10 ) = x24;
    J.insert( 24, 21 ) = -1;
    J.insert( 24, 23 ) = x11;
    J.insert( 24, 24 ) = -1;

    J.insert( 25, 10 ) = 2 * ( x11 * x12 * x14 ) * x12 * x14;
    J.insert( 25, 11 ) = 2 * ( x11 * x12 * x14 ) * x11 * x14;
    J.insert( 25, 13 ) = 2 * ( x11 * x12 * x14 ) * x11 * x12;
    J.insert( 25, 25 ) = -1;

    J.insert( 26, 0 )  = x5 / 10000;
    J.insert( 26, 4 )  = x1 / 10000;
    J.insert( 26, 26 ) = -1;

    J.insert( 27, 5 )  = 2 * ( x6 - x8 ) / x22;
    J.insert( 27, 7 )  = -2 * ( x6 - x8 ) / x22;
    J.insert( 27, 21 ) = -power2( ( x6 - x8 ) / x22 );
    J.insert( 27, 27 ) = -1;

    J.insert( 28, 15 ) = x19 / x21;
    J.insert( 28, 18 ) = x16 / x21;
    J.insert( 28, 20 ) = -x16 * x19 / power2( x21 );
    J.insert( 28, 28 ) = -1;

    J.insert( 29, 21 ) = 1;
    J.insert( 29, 26 ) = x28;
    J.insert( 29, 27 ) = x27 + 0.5 / sqrt( x28 );
    J.insert( 29, 29 ) = -1;

    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 1.5, 1, 50, 50, 500, 5000, 0, 5000, 0, 5000, 0.5, 50, 5, 0.5, 50, 50, 5, 50, 50, 5, 50, 500, 5, 50, -50, 5000,
      0.05, 50000, 5, 5000;
  }

  virtual void check_if_admissible( Vector const & x ) const override
  {
    real_type x1  = x( 0 );
    real_type x4  = x( 3 );
    real_type x5  = x( 4 );
    real_type x6  = x( 5 );
    real_type x8  = x( 7 );
    real_type x10 = x( 9 );
    real_type x22 = x( 21 );
    real_type x28 = x( 27 );
    UTILS_ASSERT( x1 > 0, "x1" );
    UTILS_ASSERT( x4 > 0, "x4" );
    UTILS_ASSERT( x5 > 0, "x5" );
    UTILS_ASSERT( x6 > 0, "x6" );
    UTILS_ASSERT( x8 > 0, "x8" );
    UTILS_ASSERT( x10 > 0, "x10" );
    UTILS_ASSERT( x22 > 0, "x22" );
    UTILS_ASSERT( x28 > 0, "x28" );
  }

  virtual void bounding_box( Vector & L, Vector & U ) const override
  {
    U.fill( real_max );
    L.fill( -real_max );
    L[0] = L[3] = L[4] = L[5] = L[7] = L[9] = L[21] = L[27] = 0;
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/
