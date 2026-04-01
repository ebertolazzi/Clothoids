/*\
 |
 |  Author:
 |    Enrico Bertolazzi
 |    University of Trento
 |    Department of Industrial Engineering
 |    Via Sommarive 9, I-38123, Povo, Trento, Italy
 |    email: enrico.bertolazzi@unitn.it
\*/

#define DENNIS_AND_GAY_BIBTEX                                                \
  "@techreport{Dennis:1983,\n"                                               \
  "  author = {J. E. Dennis, D. M. Gay, and P. A. Vu},\n"                    \
  "  title  = {A new nonlinear equations test problem},\n"                   \
  "  number = {Technical Report 83-16, Mathematical Sciences Department},\n" \
  "  year   = {1983}\n"                                                      \
  "  note   = {Rice University (1983 - revised 1985)}\n"                     \
  "}\n"

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class test_dennis_gay_6 : public NonlinearSystem
{
  real_type const summx;
  real_type const summy;
  real_type const suma;
  real_type const sumb;
  real_type const sumc;
  real_type const sumd;
  real_type const sume;
  real_type const sumf;
  real_type const x0;
  real_type const x1;
  real_type const x2;
  real_type const x3;
  real_type const x4;
  real_type const x5;

public:
  test_dennis_gay_6(
    string const & n,
    real_type      _summx,
    real_type      _summy,
    real_type      _suma,
    real_type      _sumb,
    real_type      _sumc,
    real_type      _sumd,
    real_type      _sume,
    real_type      _sumf,
    real_type      _x0,
    real_type      _x1,
    real_type      _x2,
    real_type      _x3,
    real_type      _x4,
    real_type      _x5 )
    : NonlinearSystem( n, DENNIS_AND_GAY_BIBTEX, 6 )
    , summx( _summx )
    , summy( _summy )
    , suma( _suma )
    , sumb( _sumb )
    , sumc( _sumc )
    , sumd( _sumd )
    , sume( _sume )
    , sumf( _sumf )
    , x0( _x0 )
    , x1( _x1 )
    , x2( _x2 )
    , x3( _x3 )
    , x4( _x4 )
    , x5( _x5 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type a = x( 0 );
    real_type b = summx - a;
    real_type c = x( 1 );
    real_type d = summy - c;
    real_type t = x( 2 );
    real_type u = x( 3 );
    real_type v = x( 4 );
    real_type w = x( 5 );

    // real_type tv    = t*v;
    real_type tt    = t * t;
    real_type vv    = v * v;
    real_type tsvs  = tt - vv;
    real_type ts3vs = tt - 3 * vv;
    real_type vs3ts = vv - 3 * tt;
    // real_type uw    = u*w;
    real_type uu    = u * u;
    real_type ww    = w * w;
    real_type usws  = uu - ww;
    real_type us3ws = uu - 3 * ww;
    real_type ws3us = ww - 3 * uu;

    f( 0 ) = t * a + u * b - v * c - w * d - suma;
    f( 1 ) = v * a + w * b + t * c + u * d - sumb;
    f( 2 ) = a * tsvs - 2 * c * t * v + b * usws - 2 * d * u * w - sumc;
    f( 3 ) = c * tsvs + 2 * a * t * v + d * usws + 2 * b * u * w - sumd;
    f( 4 ) = a * t * ts3vs + c * v * vs3ts + b * u * us3ws + d * w * ws3us - sume;
    f( 5 ) = c * t * ts3vs - a * v * vs3ts + d * u * us3ws - b * w * ws3us - sumf;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    real_type a = x( 0 );
    real_type b = summx - a;
    real_type c = x( 1 );
    real_type d = summy - c;
    real_type t = x( 2 );
    real_type u = x( 3 );
    real_type v = x( 4 );
    real_type w = x( 5 );

    real_type tv    = t * v;
    real_type tt    = t * t;
    real_type vv    = v * v;
    real_type tsvs  = tt - vv;
    real_type ts3vs = tt - 3 * vv;
    real_type vs3ts = vv - 3 * tt;
    real_type uw    = u * w;
    real_type uu    = u * u;
    real_type ww    = w * w;
    real_type usws  = uu - ww;
    real_type us3ws = uu - 3 * ww;
    real_type ws3us = ww - 3 * uu;

    J.insert( 0, 0 ) = t - u;
    J.insert( 0, 1 ) = -v + w;
    J.insert( 0, 2 ) = a;
    J.insert( 0, 3 ) = b;
    J.insert( 0, 4 ) = -c;
    J.insert( 0, 5 ) = -d;

    J.insert( 1, 0 ) = v - w;
    J.insert( 1, 1 ) = t - u;
    J.insert( 1, 2 ) = c;
    J.insert( 1, 3 ) = d;
    J.insert( 1, 4 ) = a;
    J.insert( 1, 5 ) = b;

    J.insert( 2, 0 ) = tsvs - usws;
    J.insert( 2, 1 ) = -2 * ( tv - uw );
    J.insert( 2, 2 ) = 2 * ( a * t - c * v );
    J.insert( 2, 3 ) = 2 * ( b * u - d * w );
    J.insert( 2, 4 ) = -2 * ( a * v + c * t );
    J.insert( 2, 5 ) = -2 * ( b * w + d * u );

    J.insert( 3, 0 ) = 2 * ( tv - uw );
    J.insert( 3, 1 ) = tsvs - usws;
    J.insert( 3, 2 ) = 2 * ( c * t + a * v );
    J.insert( 3, 3 ) = 2 * ( d * u + b * w );
    J.insert( 3, 4 ) = 2 * ( a * t - c * v );
    J.insert( 3, 5 ) = 2 * ( b * u - d * w );

    J.insert( 4, 0 ) = t * ts3vs - u * us3ws;
    J.insert( 4, 1 ) = v * vs3ts - w * ws3us;
    J.insert( 4, 2 ) = 3 * ( a * tsvs - 2 * c * tv );
    J.insert( 4, 3 ) = 3 * ( b * usws - 2 * d * uw );
    J.insert( 4, 4 ) = -3 * ( c * tsvs + 2 * a * tv );
    J.insert( 4, 5 ) = -3 * ( d * usws + 2 * b * uw );

    J.insert( 5, 0 ) = -v * vs3ts + w * ws3us;
    J.insert( 5, 1 ) = t * ts3vs - u * us3ws;
    J.insert( 5, 2 ) = 3 * ( c * tsvs + 2 * a * tv );
    J.insert( 5, 3 ) = 3 * ( d * usws + 2 * b * uw );
    J.insert( 5, 4 ) = 3 * ( a * tsvs - 2 * c * tv );
    J.insert( 5, 5 ) = 3 * ( b * usws - 2 * d * uw );

    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    x_vec[0].resize( n );
    x_vec[0] << x0, x1, x2, x3, x4, x5;
  }

  virtual void check_if_admissible( Vector const & x ) const override
  {
    for ( integer i = 0; i < n; ++i ) UTILS_ASSERT( std::abs( x( i ) ) < 100, "Bad range" );
  }

  virtual void bounding_box( Vector & L, Vector & U ) const override
  {
    U.fill( 100 );
    L.fill( -100 );
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class DennisAndGay6eqN1 : public test_dennis_gay_6
{
public:
  DennisAndGay6eqN1()
    : test_dennis_gay_6(
        "Dennis and Gay 6 eq N 1",
        0.485,
        -0.0019,
        -0.0581,
        0.015,
        0.105,
        0.0406,
        0.167,
        -0.399,
        0.299,
        -0.0273,
        -0.474,
        0.474,
        -0.0892,
        0.0892 )
  {
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class DennisAndGay6eqN2 : public test_dennis_gay_6
{
public:
  DennisAndGay6eqN2()
    : test_dennis_gay_6(
        "Dennis and Gay 6 eq N 2",
        -0.69,
        -0.044,
        -1.57,
        -1.31,
        -2.65,
        2.0,
        -12.6,
        9.48,
        -0.3,
        0.3,
        -1.2,
        2.69,
        1.59,
        -1.5 )
  {
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class DennisAndGay6eqN3 : public test_dennis_gay_6
{
public:
  DennisAndGay6eqN3()
    : test_dennis_gay_6(
        "Dennis and Gay 6 eq N 3",
        -0.816,
        -0.017,
        -1.826,
        -0.754,
        -4.839,
        -3.259,
        -14.023,
        15.467,
        -0.041,
        0.03,
        -2.565,
        2.565,
        -0.754,
        0.754 )
  {
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class DennisAndGay6eqN4 : public test_dennis_gay_6
{
public:
  DennisAndGay6eqN4()
    : test_dennis_gay_6(
        "Dennis and Gay 6 eq N 4",
        -0.809,
        -0.021,
        -2.04,
        -0.614,
        -6.903,
        -2.934,
        -26.328,
        18.639,
        -0.056,
        0.026,
        -2.991,
        2.991,
        -0.568,
        0.568 )
  {
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class DennisAndGay6eqN5 : public test_dennis_gay_6
{
public:
  DennisAndGay6eqN5()
    : test_dennis_gay_6(
        "Dennis and Gay 6 eq N 5",
        -0.807,
        -0.021,
        -2.379,
        -0.364,
        -10.541,
        -1.961,
        -51.551,
        21.053,
        -0.074,
        0.013,
        -3.632,
        3.632,
        -0.289,
        0.289 )
  {
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class test_dennis_gay_8 : public NonlinearSystem
{
  real_type const summx;
  real_type const summy;
  real_type const suma;
  real_type const sumb;
  real_type const sumc;
  real_type const sumd;
  real_type const sume;
  real_type const sumf;
  real_type const x0;
  real_type const x1;
  real_type const x2;
  real_type const x3;
  real_type const x4;
  real_type const x5;
  real_type const x6;
  real_type const x7;

public:
  test_dennis_gay_8(
    string const & n,
    real_type      _summx,
    real_type      _summy,
    real_type      _suma,
    real_type      _sumb,
    real_type      _sumc,
    real_type      _sumd,
    real_type      _sume,
    real_type      _sumf,
    real_type      _x0,
    real_type      _x1,
    real_type      _x2,
    real_type      _x3,
    real_type      _x4,
    real_type      _x5,
    real_type      _x6,
    real_type      _x7 )
    : NonlinearSystem( n, DENNIS_AND_GAY_BIBTEX, 8 )
    , summx( _summx )
    , summy( _summy )
    , suma( _suma )
    , sumb( _sumb )
    , sumc( _sumc )
    , sumd( _sumd )
    , sume( _sume )
    , sumf( _sumf )
    , x0( _x0 )
    , x1( _x1 )
    , x2( _x2 )
    , x3( _x3 )
    , x4( _x4 )
    , x5( _x5 )
    , x6( _x6 )
    , x7( _x7 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type a = x( 0 );
    real_type b = x( 1 );
    real_type c = x( 2 );
    real_type d = x( 3 );
    real_type t = x( 4 );
    real_type u = x( 5 );
    real_type v = x( 6 );
    real_type w = x( 7 );

    // real_type tv    = t*v;
    real_type tt    = t * t;
    real_type vv    = v * v;
    real_type tsvs  = tt - vv;
    real_type ts3vs = tt - 3 * vv;
    real_type vs3ts = vv - 3 * tt;
    // real_type uw    = u*w;
    real_type uu    = u * u;
    real_type ww    = w * w;
    real_type usws  = uu - ww;
    real_type us3ws = uu - 3 * ww;
    real_type ws3us = ww - 3 * uu;

    f( 0 ) = a + b - summx;
    f( 1 ) = c + d - summy;
    f( 2 ) = t * a + u * b - v * c - w * d - suma;
    f( 3 ) = v * a + w * b + t * c + u * d - sumb;
    f( 4 ) = a * tsvs - 2 * c * t * v + b * usws - 2 * d * u * w - sumc;
    f( 5 ) = c * tsvs + 2 * a * t * v + d * usws + 2 * b * u * w - sumd;
    f( 6 ) = a * t * ts3vs + c * v * vs3ts + b * u * us3ws + d * w * ws3us - sume;
    f( 7 ) = c * t * ts3vs - a * v * vs3ts + d * u * us3ws - b * w * ws3us - sumf;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    real_type a = x( 0 );
    real_type b = x( 1 );
    real_type c = x( 2 );
    real_type d = x( 3 );
    real_type t = x( 4 );
    real_type u = x( 5 );
    real_type v = x( 6 );
    real_type w = x( 7 );

    real_type tv    = t * v;
    real_type tt    = t * t;
    real_type vv    = v * v;
    real_type tsvs  = tt - vv;
    real_type ts3vs = tt - 3 * vv;
    real_type vs3ts = vv - 3 * tt;
    real_type uw    = u * w;
    real_type uu    = u * u;
    real_type ww    = w * w;
    real_type usws  = uu - ww;
    real_type us3ws = uu - 3 * ww;
    real_type ws3us = ww - 3 * uu;

    J.insert( 0, 0 ) = 1;
    J.insert( 0, 1 ) = 1;
    J.insert( 0, 2 ) = 0;
    J.insert( 0, 3 ) = 0;
    J.insert( 0, 4 ) = 0;
    J.insert( 0, 5 ) = 0;
    J.insert( 0, 6 ) = 0;
    J.insert( 0, 7 ) = 0;

    J.insert( 1, 0 ) = 0;
    J.insert( 1, 1 ) = 0;
    J.insert( 1, 2 ) = 1;
    J.insert( 1, 3 ) = 1;
    J.insert( 1, 4 ) = 0;
    J.insert( 1, 5 ) = 0;
    J.insert( 1, 6 ) = 0;
    J.insert( 1, 7 ) = 0;

    J.insert( 2, 0 ) = t;
    J.insert( 2, 1 ) = u;
    J.insert( 2, 2 ) = -v;
    J.insert( 2, 3 ) = -w;
    J.insert( 2, 4 ) = a;
    J.insert( 2, 5 ) = b;
    J.insert( 2, 6 ) = -c;
    J.insert( 2, 7 ) = -d;

    J.insert( 3, 0 ) = v;
    J.insert( 3, 1 ) = w;
    J.insert( 3, 2 ) = t;
    J.insert( 3, 3 ) = u;
    J.insert( 3, 4 ) = c;
    J.insert( 3, 5 ) = d;
    J.insert( 3, 6 ) = a;
    J.insert( 3, 7 ) = b;

    J.insert( 4, 0 ) = tsvs;
    J.insert( 4, 1 ) = usws;
    J.insert( 4, 2 ) = -2 * tv;
    J.insert( 4, 3 ) = -2 * uw;
    J.insert( 4, 4 ) = 2 * ( a * t - c * v );
    J.insert( 4, 5 ) = 2 * ( b * u - d * w );
    J.insert( 4, 6 ) = -2 * ( a * v + c * t );
    J.insert( 4, 7 ) = -2 * ( b * w + d * u );

    J.insert( 5, 0 ) = 2 * tv;
    J.insert( 5, 1 ) = 2 * uw;
    J.insert( 5, 2 ) = tsvs;
    J.insert( 5, 3 ) = usws;
    J.insert( 5, 4 ) = 2 * ( c * t + a * v );
    J.insert( 5, 5 ) = 2 * ( d * u + b * w );
    J.insert( 5, 6 ) = 2 * ( a * t - c * v );
    J.insert( 5, 7 ) = 2 * ( b * u - d * w );

    J.insert( 6, 0 ) = t * ts3vs;
    J.insert( 6, 1 ) = u * us3ws;
    J.insert( 6, 2 ) = v * vs3ts;
    J.insert( 6, 3 ) = w * ws3us;
    J.insert( 6, 4 ) = 3 * ( a * tsvs - 2 * c * tv );
    J.insert( 6, 5 ) = 3 * ( b * usws - 2 * d * uw );
    J.insert( 6, 6 ) = -3 * ( c * tsvs + 2 * a * tv );
    J.insert( 6, 7 ) = -3 * ( d * usws + 2 * b * uw );

    J.insert( 7, 0 ) = -v * vs3ts;
    J.insert( 7, 1 ) = -w * ws3us;
    J.insert( 7, 2 ) = t * ts3vs;
    J.insert( 7, 3 ) = u * us3ws;
    J.insert( 7, 4 ) = 3 * ( c * tsvs + 2 * a * tv );
    J.insert( 7, 5 ) = 3 * ( d * usws + 2 * b * uw );
    J.insert( 7, 6 ) = 3 * ( a * tsvs - 2 * c * tv );
    J.insert( 7, 7 ) = 3 * ( b * usws - 2 * d * uw );

    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    x_vec[0].resize( n );
    x_vec[0] << x0, x1, x2, x3, x4, x5, x6, x7;
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class DennisAndGay8eqN1 : public test_dennis_gay_8
{
public:
  DennisAndGay8eqN1()
    : test_dennis_gay_8(
        "Dennis and Gay 8 eq N 1",
        0.485,
        -0.0019,
        -0.0581,
        0.015,
        0.105,
        0.0406,
        0.167,
        -0.399,
        0.299,
        0.186,
        -0.0273,
        0.0254,
        -0.474,
        0.474,
        -0.0892,
        0.0892 )
  {
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class DennisAndGay8eqN2 : public test_dennis_gay_8
{
public:
  DennisAndGay8eqN2()
    : test_dennis_gay_8(
        "Dennis and Gay 8 eq N 2",
        -0.69,
        -0.044,
        -1.57,
        -1.31,
        -2.65,
        2.0,
        -12.6,
        9.48,
        -0.3,
        -0.39,
        0.3,
        -0.344,
        -1.2,
        2.69,
        1.59,
        -1.5 )
  {
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class DennisAndGay8eqN3 : public test_dennis_gay_8
{
public:
  DennisAndGay8eqN3()
    : test_dennis_gay_8(
        "Dennis and Gay 8 eq N 3",
        -0.816,
        -0.017,
        -1.826,
        -0.754,
        -4.839,
        -3.259,
        -14.023,
        15.467,
        -0.041,
        -0.775,
        0.03,
        -0.047,
        -2.565,
        2.565,
        -0.754,
        0.754 )
  {
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class DennisAndGay8eqN4 : public test_dennis_gay_8
{
public:
  DennisAndGay8eqN4()
    : test_dennis_gay_8(
        "Dennis and Gay 8 eq N 4",
        -0.809,
        -0.021,
        -2.04,
        -0.614,
        -6.903,
        -2.934,
        -26.328,
        18.639,
        -0.056,
        -0.753,
        0.026,
        -0.047,
        -2.991,
        2.991,
        -0.568,
        0.568 )
  {
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class DennisAndGay8eqN5 : public test_dennis_gay_8
{
public:
  DennisAndGay8eqN5()
    : test_dennis_gay_8(
        "Dennis and Gay 8 eq N 5",
        -0.807,
        -0.021,
        -2.379,
        -0.364,
        -10.541,
        -1.961,
        -51.551,
        21.053,
        -0.074,
        -0.733,
        0.013,
        -0.034,
        -3.632,
        3.632,
        -0.289,
        0.289 )
  {
  }
};
