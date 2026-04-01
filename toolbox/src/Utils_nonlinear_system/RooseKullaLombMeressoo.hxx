/*\
 |
 |  Author:
 |    Enrico Bertolazzi
 |    University of Trento
 |    Department of Industrial Engineering
 |    Via Sommarive 9, I-38123, Povo, Trento, Italy
 |    email: enrico.bertolazzi@unitn.it
\*/

/*
"@Article{Spedicato1997,\n"
"  author  = {Spedicato, E. and Huang, Z.},\n"
"  title   = {Numerical experience with newton-like methods\n"
"             for nonlinear algebraic systems},\n"
"  journal = {Computing},\n"
"  year    = {1997},\n"
"  volume  = {58},\n"
"  number  = {1},\n"
"  pages   = {69--89},\n"
"  doi     = {10.1007/BF02684472},\n"
"}\n"

  Numerical Experience with Newton-like Methods for Nonlinear Algebraic Systems
  E. Spedicato, Z. Huang
  COmputing N. 58, 1997

  Test Problems for Unconstrained Optimization
  Ladislav Luksan Jan Vlcek
  Technical report No. 897
  January 2003

  Sparse Test Problems for Unconstrained Optimization
  Ladislav Luksan Jan Vlcek
  Technical report No. 1064
  January 2010

  A MODIFIED NEWTON METHOD FOR SOLVING NON-LINEAR ALGEBRAIC EQUATIONS
  Satya N. Atluri*, Chein-Shan Liu**, and Chung-Lun Kuo***
  Journal of Marine Science and Technology, Vol. 17, No. 3, pp. 238-247 (2009)
 */

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

#define RKM_BIBTEX                                                           \
  "@book{meresoo:1990,\n"                                                    \
  "  title     = {Test Examples of Systems of Nonlinear Equations: Version " \
  "3-90},\n"                                                                 \
  "  author    = {Meresoo, T. and Roose, A. and Kulla,\n"                    \
  "               V. and Estonian Software and Computer Service Company},\n" \
  "  year      = 1990,\n"                                                    \
  "  publisher = {Estonian Software and Computer Service Company}\n"         \
  "}\n"

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo129 : public NonlinearSystem
{
public:
  RooseKullaLombMeressoo129() : NonlinearSystem( "Roose Kulla Lomb Meressoo N.129", RKM_BIBTEX, 3 ) {}

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = x( 0 ) + x( 1 ) - 2;
    f( 1 ) = x( 0 ) - log( x( 1 ) ) + x( 2 ) - 2;
    f( 2 ) = x( 1 ) * x( 1 ) - 2 * x( 2 ) + 1;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    J.insert( 0, 0 ) = 1;
    J.insert( 0, 1 ) = 1;
    J.insert( 0, 2 ) = 0;

    J.insert( 1, 0 ) = 1;
    J.insert( 1, 1 ) = -1 / x( 1 );
    J.insert( 1, 2 ) = 1;

    J.insert( 2, 0 ) = 0;
    J.insert( 2, 1 ) = 2 * x( 1 );
    J.insert( 2, 2 ) = -2;

    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 2 );
    auto & x0{ x_vec[0] };
    auto & x1{ x_vec[1] };
    x0.resize( n );
    x1.resize( n );
    x0.fill( 1 );
    x1 << -0.285888702509893511, 2.28588870250989329, 3.11264358013118247;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 0.5 );
  }

  string note() const { return "Jacobian is singular at the solution"; }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo201 : public NonlinearSystem
{
public:
  RooseKullaLombMeressoo201( integer neq ) : NonlinearSystem( "Roose Kulla Lomb Meressoo N.201", RKM_BIBTEX, neq )
  {
    check_min_equations( n, 1 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = 1 - x( 0 );
    for ( integer k = 1; k < n; ++k ) f( k ) = 10 * k * power2( x( k ) - x( k - 1 ) );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    J.insert( 0, 0 ) = -1;
    for ( integer k = 1; k < n; ++k )
    {
      J.insert( k, k - 1 ) = -20 * k * ( x( k ) - x( k - 1 ) );
      J.insert( k, k )     = 20 * k * ( x( k ) - x( k - 1 ) );
    }
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 1 );
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( -1.2 );
    x0( n - 1 ) = -1;
  }

  string note() const { return "Jacobian is singular at the solution"; }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo202 : public NonlinearSystem
{
public:
  RooseKullaLombMeressoo202( integer neq ) : NonlinearSystem( "Roose Kulla Lomb Meressoo N.202", RKM_BIBTEX, neq )
  {
    check_min_equations( n, 1 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    for ( integer k = 0; k < n - 1; ++k ) f( k ) = x( k ) - 0.1 * power2( x( k + 1 ) );
    f( n - 1 ) = x( n - 1 ) - 0.1 * power2( x( 0 ) );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer k = 0; k < n - 1; ++k )
    {
      J.insert( k, k )     = 1;
      J.insert( k, k + 1 ) = -0.2 * x( k + 1 );
    }
    J.insert( n - 1, n - 1 ) = 1;
    J.insert( n - 1, 0 )     = -0.2 * x( 0 );
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 2 );
    auto & x0{ x_vec[0] };
    auto & x1{ x_vec[1] };
    x0.resize( n );
    x1.resize( n );
    x0.setZero();
    x1.fill( 10 );
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 2 );
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo203 : public NonlinearSystem
{
  using Matrix = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;

public:
  RooseKullaLombMeressoo203( integer neq ) : NonlinearSystem( "Roose Kulla Lomb Meressoo N.203", RKM_BIBTEX, neq )
  {
    check_min_equations( n, 2 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type sum = 0;
    for ( integer i = 0; i < n; ++i ) sum += x( i );
    for ( integer k = 0; k < n - 1; ++k ) f( k ) = x( k ) - ( n + 1 ) + sum;
    f( n - 1 ) = sum - 1;
  }

  virtual void jacobian( Vector const &, SparseMatrix & J ) const override
  {
    Matrix J_full( n, n );
    J_full.fill( 1 );
    for ( integer i = 0; i < n - 1; ++i ) J_full( i, i ) += 1;
    J.resize( n, n );
    J = J_full.sparseView();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    switch ( n )
    {
      case 2:
      {
        x_vec.resize( 2 );
        auto & x0{ x_vec[0] };
        auto & x1{ x_vec[1] };
        x0.resize( n );
        x1.resize( n );
        x0.fill( 1 );
        x1 << 0.5, 2;
      }
      break;
      case 5:
      {
        x_vec.resize( 3 );
        auto & x0{ x_vec[0] };
        auto & x1{ x_vec[1] };
        auto & x2{ x_vec[2] };
        x0.resize( n );
        x1.resize( n );
        x2.resize( n );
        x0.fill( 1 );
        x1.fill( 0.91635458253384926 );
        x1( 4 ) = 1.4182270873307543;
        x2.fill( -5.7904308849411582E-1 );
        x2( 4 ) = 9.8952154424705789;
      }
      break;
      case 10:
      case 20:
      case 30:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0.fill( 1 );
      }
    }
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 0.5 );
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo204 : public NonlinearSystem
{
public:
  RooseKullaLombMeressoo204( integer neq ) : NonlinearSystem( "Roose Kulla Lomb Meressoo N.204", RKM_BIBTEX, neq )
  {
    check_even( n, 1 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    for ( integer k = 0; k < n; k += 2 ) f( k ) = 1 - x( k );
    for ( integer k = 1; k < n; k += 2 ) f( k ) = 10 * ( x( k ) - power2( x( k - 1 ) ) );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer k = 0; k < n; k += 2 ) J.insert( k, k ) = -1;
    for ( integer k = 1; k < n; k += 2 )
    {
      J.insert( k, k - 1 ) = -20 * x( k - 1 );
      J.insert( k, k )     = 10;
    }
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 1 );
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    for ( integer k{ 0 }; k < n; k += 2 ) x0( k ) = -1.2;
    for ( integer k{ 1 }; k < n; k += 2 ) x0( k ) = 1;
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo205 : public NonlinearSystem
{
public:
  RooseKullaLombMeressoo205( integer neq ) : NonlinearSystem( "Roose Kulla Lomb Meressoo N.205", RKM_BIBTEX, neq )
  {
    check_min_equations( n, 1 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type acc = 0;
    for ( integer j = 0; j < n; ++j ) acc += power3( x( j ) );
    acc /= 2 * n;
    for ( integer j = 0; j < n; ++j ) f( j ) = x( j ) - acc - ( 0.5 / n ) * ( j + 1 );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type bf = -1.5 / n;
    for ( integer i = 0; i < n; ++i )
    {
      for ( integer j = 0; j < n; ++j )
      {
        real_type tmp = bf * power2( x( j ) );
        if ( i == j ) tmp += 1;
        J.insert( i, j ) = tmp;
      }
    }
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    switch ( n )
    {
      case 2:
      {
        x_vec.resize( 3 );
        auto & x0{ x_vec[0] };
        auto & x1{ x_vec[1] };
        auto & x2{ x_vec[2] };
        x0.resize( n );
        x1.resize( n );
        x2.resize( n );
        x0 << 1.01242808774686965, 1.26242808774686965;
        x1 << -1.68508582524678752, -1.43508582524678752;
        x2 << 0.297657737499917863, 0.547657737499917863;
      }
      break;

      case 5:
      {
        x_vec.resize( 3 );
        auto & x0{ x_vec[0] };
        auto & x1{ x_vec[1] };
        auto & x2{ x_vec[2] };
        x0.resize( n );
        x1.resize( n );
        x2.resize( n );
        x0 << 1, 1.1, 1.2, 1.3, 1.4;
        x1 << -1.72736184954957039, -1.62736184954957053, -1.52736184954957044, -1.42736184954957057,
          -1.32736184954957048;
        x2 << 0.127361849549570388, 0.227361849549570394, 0.327361849549570427, 0.427361849549570405,
          0.527361849549570438;
      }
      break;

      case 10:
      {
        x_vec.resize( 3 );
        auto & x0{ x_vec[0] };
        auto & x1{ x_vec[1] };
        auto & x2{ x_vec[2] };
        x0.resize( n );
        x1.resize( n );
        x2.resize( n );
        x0 << 0.994471118343573157, 1.04447111834357309, 1.09447111834357313, 1.14447111834357318, 1.194471118343573,
          1.24447111834357305, 1.29447111834357309, 1.34447111834357313, 1.39447111834357318, 1.444471118343573;
        x1 << -1.74181474184792218, -1.69181474184792235, -1.64181474184792231, -1.59181474184792227,
          -1.54181474184792222, -1.49181474184792218, -1.44181474184792213, -1.39181474184792231, -1.34181474184792227,
          -1.29181474184792222;
        x2 << 0.0723436235043492942, 0.122343623504349297, 0.172343623504349314, 0.222343623504349303,
          0.272343623504349264, 0.322343623504349364, 0.372343623504349297, 0.422343623504349341, 0.472343623504349275,
          0.522343623504349264;
      }
      break;

      case 20:
      {
        x_vec.resize( 3 );
        auto & x0{ x_vec[0] };
        auto & x1{ x_vec[1] };
        auto & x2{ x_vec[2] };
        x0.resize( n );
        x1.resize( n );
        x2.resize( n );
        x0 << 0.991518329652893882, 1.0165183296528939, 1.04151832965289382, 1.06651832965289395, 1.09151832965289386,
          1.11651832965289377, 1.1415183296528939, 1.16651832965289382, 1.19151832965289395, 1.21651832965289386,
          1.24151832965289377, 1.2665183296528939, 1.29151832965289382, 1.31651832965289395, 1.34151832965289386,
          1.36651832965289377, 1.3915183296528939, 1.41651832965289382, 1.44151832965289395, 1.46651832965289386;

        x1 << -1.74911100354538962, -1.72411100354538971, -1.69911100354538958, -1.67411100354538966,
          -1.64911100354538975, -1.62411100354538962, -1.59911100354538971, -1.57411100354538958, -1.54911100354538966,
          -1.52411100354538975, -1.49911100354538962, -1.47411100354538971, -1.44911100354538958, -1.42411100354538966,
          -1.39911100354538975, -1.37411100354538962, -1.34911100354538971, -1.32411100354538958, -1.29911100354538966,
          -1.27411100354538975;
        x2 << 0.0450926738924959589, 0.0700926738924959603, 0.0950926738924959686, 0.120092673892495963,
          0.145092673892495971, 0.170092673892495994, 0.19509267389249596, 0.220092673892495982, 0.245092673892495949,
          0.270092673892495971, 0.295092673892495994, 0.320092673892496016, 0.345092673892495982, 0.370092673892496005,
          0.395092673892495971, 0.420092673892495994, 0.445092673892496016, 0.470092673892495982, 0.495092673892496005,
          0.520092673892496027;
      }
      break;

      case 30:
      {
        x_vec.resize( 3 );
        auto & x0{ x_vec[0] };
        auto & x1{ x_vec[1] };
        auto & x2{ x_vec[2] };
        x0.resize( n );
        x1.resize( n );
        x2.resize( n );

        x0 << 0.990509030546815938, 1.00717569721348266, 1.02384236388014926, 1.04050903054681587, 1.0571756972134827,
          1.07384236388014931, 1.09050903054681592, 1.10717569721348252, 1.12384236388014935, 1.14050903054681596,
          1.15717569721348257, 1.1738423638801494, 1.19050903054681601, 1.20717569721348261, 1.22384236388014922,
          1.24050903054681605, 1.25717569721348266, 1.27384236388014926, 1.29050903054681587, 1.3071756972134827,
          1.32384236388014931, 1.34050903054681592, 1.35717569721348252, 1.37384236388014935, 1.39050903054681596,
          1.40717569721348257, 1.4238423638801494, 1.44050903054681601, 1.45717569721348261, 1.47384236388014922;

        x1 << -1.75155356027278342, -1.73488689360611681, -1.71822022693944998, -1.70155356027278337,
          -1.68488689360611676, -1.66822022693945016, -1.65155356027278333, -1.63488689360611672, -1.61822022693945011,
          -1.60155356027278351, -1.58488689360611668, -1.56822022693945007, -1.55155356027278346, -1.53488689360611685,
          -1.51822022693945002, -1.50155356027278342, -1.48488689360611681, -1.4682202269394502, -1.45155356027278337,
          -1.43488689360611676, -1.41822022693945016, -1.40155356027278355, -1.38488689360611672, -1.36822022693945011,
          -1.35155356027278351, -1.33488689360611668, -1.31822022693945007, -1.30155356027278346, -1.28488689360611685,
          -1.26822022693945002;

        x2 << 0.0360445297259673891, 0.0527111963926340521, 0.0693778630593007289, 0.0860445297259673919,
          0.102711196392634055, 0.119377863059300732, 0.136044529725967395, 0.152711196392634058, 0.169377863059300721,
          0.186044529725967384, 0.202711196392634047, 0.219377863059300737, 0.2360445297259674, 0.252711196392634063,
          0.269377863059300726, 0.286044529725967389, 0.302711196392634052, 0.319377863059300715, 0.336044529725967378,
          0.352711196392634041, 0.369377863059300704, 0.386044529725967367, 0.40271119639263403, 0.419377863059300748,
          0.436044529725967411, 0.452711196392634074, 0.469377863059300737, 0.4860445297259674, 0.502711196392634063,
          0.519377863059300782;
      }
      break;
    }
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 1.5 );
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo206 : public NonlinearSystem
{
  real_type h2;

public:
  RooseKullaLombMeressoo206( integer neq ) : NonlinearSystem( "Roose Kulla Lomb Meressoo N.206", RKM_BIBTEX, neq )
  {
    check_min_equations( n, 1 );
    h2 = 1.0 / power2( n + 1 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    for ( integer k = 0; k < n; ++k )
    {
      real_type xc = x( k );
      real_type xp = k < n - 1 ? x( k + 1 ) : 0;
      real_type xm = k > 0 ? x( k - 1 ) : 0;
      f( k )       = xp - 2 * xc + xm - h2 * exp( xc );
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer k = 1; k < n; ++k ) J.insert( k, k - 1 ) = 1;
    for ( integer k = 0; k < n - 1; ++k ) J.insert( k, k + 1 ) = 1;
    for ( integer k = 0; k < n; ++k ) J.insert( k, k ) = -2 - h2 * exp( x( k ) );
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    switch ( n )
    {
      case 2:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << -0.100488400337317069, -0.100488400337317069;
      }
      break;
      case 5:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << -0.0635730237796020142, -0.101079225590438845, -0.113478165704209086, -0.101079225590438845,
          -0.0635730237796020142;
      }
      break;
      case 10:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << -0.0380470822832765579, -0.068138233862133038, -0.0905092917544646769, -0.10533104513306768,
          -0.112714562777264757, -0.112714562777264757, -0.10533104513306768, -0.0905092917544646769,
          -0.068138233862133038, -0.038047082283276551;
      }
      break;
      case 20:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << -0.0209484001798467684, -0.0396762346150324255, -0.0562247027012899511, -0.0706295728406963474,
          -0.0829215019487556798, -0.0931263031751204673, -0.101265166939055559, -0.107354839384380743,
          -0.11140776149579644, -0.113432171358331196, -0.113432171358331182, -0.111407761495796426,
          -0.107354839384380729, -0.101265166939055545, -0.0931263031751204534, -0.0829215019487556521,
          -0.0706295728406963197, -0.0562247027012899234, -0.0396762346150324047, -0.020948400179846758;
      }
      break;
      case 30:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << -0.0144369806182624797, -0.0278482934603717126, -0.0402476022487183935, -0.0516473776603994031,
          -0.0620589494503979225, -0.0714925529960212941, -0.0799573706823796804, -0.087461568494601391,
          -0.094012328134088427, -0.0996158749325536247, -0.104277501798100969, -0.108001589391558919,
          -0.110791622698076339, -0.112650204128129788, -0.113579063253107349, -0.113579063253107335,
          -0.112650204128129788, -0.110791622698076339, -0.108001589391558905, -0.104277501798100955,
          -0.0996158749325536108, -0.0940123281340884132, -0.0874615684946013633, -0.0799573706823796526,
          -0.0714925529960212663, -0.0620589494503978947, -0.0516473776603993753, -0.0402476022487183727,
          -0.0278482934603716953, -0.014436980618262471;
      }
      break;
    }
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.setZero();
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo207 : public NonlinearSystem
{
  real_type gamma;

public:
  RooseKullaLombMeressoo207( integer neq )
    : NonlinearSystem( "Roose Kulla Lomb Meressoo N.207", RKM_BIBTEX, neq ), gamma( 0.1 )
  {
    check_min_equations( n, 1 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    for ( integer k = 0; k < n; ++k )
    {
      real_type xc = x( k );
      real_type xp = k < n - 1 ? x( k + 1 ) : 0;
      real_type xm = k > 0 ? x( k - 1 ) : 0;
      f( k )       = ( 3 - gamma * xc ) * xc + 1 - xm - 2 * xp;
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer k = 1; k < n; ++k ) J.insert( k, k - 1 ) = -1;
    for ( integer k = 0; k < n - 1; ++k ) J.insert( k, k + 1 ) = -2;
    for ( integer k = 0; k < n; ++k ) J.insert( k, k ) = 3 - 2 * gamma * x( k );
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    switch ( n )
    {
      case 2:
      {
        x_vec.resize( 2 );
        auto & x0{ x_vec[0] };
        auto & x1{ x_vec[1] };
        x0.resize( n );
        x1.resize( n );
        x0 << 20.9686582909684525, 9.96875591028267927;
        x1 << -0.685453889794504945, -0.551673186443478292;
      }
      break;
      case 5:
      {
        x_vec.resize( 2 );
        auto & x0{ x_vec[0] };
        auto & x1{ x_vec[1] };
        x0.resize( n );
        x1.resize( n );
        x0 << 10.8597777768098975, 10.8929279971301387, 5.47670908975717818, 1.76888251337925362, 0.25852195788334581;
        x1 << -1.5293511879989905, -1.91097253481018181, -1.78437400965571991, -1.38027427739523056,
          -0.773482265306932204;
      }
      break;
      case 10:
      {
        x_vec.resize( 2 );
        auto & x0{ x_vec[0] };
        auto & x1{ x_vec[1] };
        x0.resize( n );
        x1.resize( n );
        x0 << 14.4359473997111909, 11.7340922332053701, 3.99871862308393045, -0.16845571330670267, -1.75346174786930264,
          -2.19969617021264563, -2.16474654344672679, -1.88157810993199437, -1.41701070236339088, -0.785122965109708248;
        x1 << -1.96909750815374274, -2.64751351206147811, -2.83718790384275144, -2.8345068598189691,
          -2.73488779472511823, -2.55905882466501389, -2.29858344303975581, -1.93252004445795178, -1.43622003127863795,
          -0.791206423601281572;
      }
      break;
      case 20:
      {
        x_vec.resize( 2 );
        auto & x0{ x_vec[0] };
        auto & x1{ x_vec[1] };
        x0.resize( n );
        x1.resize( n );
        x0 << 15.2764638454963269, 11.7461783871066672, 3.58240032282524856, -0.641168312964383369,
          -2.27350747113667895, -2.84811886128854175, -3.04101360876585014, -3.09984917093945933, -3.1097201961549672,
          -3.0981736936814257, -3.07233445425563545, -3.03137678448305348, -2.97036021007190376, -2.88100391174524884,
          -2.75133493955649211, -2.56499265094333007, -2.3007808716064142, -1.93335461289545485, -1.43653448650018212,
          -0.79130598984776257;
        x1 << -2.04511475622648753, -2.77679685264649745, -3.02816793889987679, -3.11234353533562258,
          -3.13876544765076693, -3.14456883057665504, -3.14188617855131014, -3.13411729048722565, -3.1213674059817329,
          -3.10213918788524357, -3.07368845588766515, -3.03184112608173084, -2.97052049186879019, -2.88105977439294048,
          -2.75135468683627105, -2.56499977369673227, -2.3007835090801767, -1.93335561455466376, -1.43653486390840945,
          -0.79130610934649992;
      }
      break;
      case 30:
      {
        x_vec.resize( 2 );
        auto & x0{ x_vec[0] };
        auto & x1{ x_vec[1] };
        x0.resize( n );
        x1.resize( n );

        x0 << 15.2931719204216012, 11.745702511253814, 3.57389143252621277, -0.650649105411801876, -2.28408658729946801,
          -2.86165790515736784, -3.05989786239379979, -3.12716658742612319, -3.14975949321800242, -3.15710518936978879,
          -3.15914369628296443, -3.15917239442777253, -3.15820525438590227, -3.15643470580650387, -3.1537034341177943,
          -3.14963006579175797, -3.14360185919570911, -3.13469938835462081, -3.12156516570161857, -3.10220650856120583,
          -3.07371144107897365, -3.0318490084888472, -2.97052321270752451, -2.88106072267857449, -2.75135502205215898,
          -2.56499989460753319, -2.30078355385205269, -1.93335563155811641, -1.43653487031502247, -0.791306111375025267;

        x1 << -2.04665682215098732, -2.77942544060933994, -3.03207003883383752, -3.11806477396577719,
          -3.14717853826406069, -3.15697205799869751, -3.1601924466152429, -3.16114345590571721, -3.16126035799173088,
          -3.16099716158523725, -3.16046071615948687, -3.15961809036597874, -3.15835610131762179, -3.15648576992995533,
          -3.15372072502463707, -3.14963592314447371, -3.14360384462249964, -3.13470006195778073, -3.12156539454732718,
          -3.10220658646387104, -3.07371146767713377, -3.03184901761026149, -2.97052321585604062, -2.88106072377591715,
          -2.75135502244006558, -2.56499989474744927, -2.30078355390386191, -1.93335563157779244, -1.4365348703224361,
          -0.791306111377372723;
      }
      break;
    }
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( -1.0 );
  }

  virtual void check_if_admissible( Vector const & x ) const override
  {
    for ( integer i = 3; i < n; ++i ) UTILS_ASSERT( x( i ) < 0, "Bad range" );
  }

  virtual void bounding_box( Vector & L, Vector & U ) const override
  {
    integer i = 0;
    for ( ; i < 3; ++i )
    {
      U[i] = real_max;
      L[i] = -real_max;
    }
    for ( ; i < n; ++i )
    {
      U[i] = 0;
      L[i] = -real_max;
    }
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo208 : public NonlinearSystem
{
  real_type const K1;
  real_type const K2;
  real_type const K3;
  integer const   R1;
  integer const   R2;

public:
  RooseKullaLombMeressoo208( integer neq )
    : NonlinearSystem( "Roose Kulla Lomb Meressoo N.208", RKM_BIBTEX, neq ), K1( 1 ), K2( 1 ), K3( 1 ), R1( 3 ), R2( 3 )
  {
    check_min_equations( n, 1 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    integer zero = 0;
    for ( integer k = 0; k < n; ++k )
    {
      real_type sum = 0;
      for ( integer i = max( zero, k - R1 ); i <= min( n - 1, k + R2 ); ++i )
        if ( i != k ) sum += x( i ) * ( 1 + x( i ) );
      f( k ) = ( K1 + K2 * x( k ) * x( k ) ) * x( k ) + 1 - K3 * sum;
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    integer zero = 0;
    for ( integer i = 0; i < n; ++i )
    {
      for ( integer j = max( zero, i - R1 ); j <= min( n - 1, i + R2 ); ++j )
      {
        if ( i == j )
          J.insert( i, j ) = K1 + 3 * K2 * x( i ) * x( i );
        else
          J.insert( i, j ) = -K3 * ( 1 + 2 * x( j ) );
      }
    }
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    switch ( n )
    {
      case 2:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << -0.754877666246692725, -0.754877666246692725;
      }
      break;
      case 5:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << -0.808025022386814507, -0.87166824608729554, -0.87166824608729554, -0.87166824608729554,
          -0.808025022386814507;
      }
      break;
      case 10:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << -0.801937946359539633, -0.842558051527937502, -0.881216883188142397, -0.911897857643132359,
          -0.891348792903065235, -0.891348792903065235, -0.911897857643132359, -0.881216883188142397,
          -0.842558051527937502, -0.801937946359539633;
      }
      break;
      case 20:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << -0.800413820399077491, -0.839125061489862767, -0.880560005368091314, -0.920665554037870759,
          -0.897476024922029492, -0.88270771158442396, -0.880029836830879542, -0.889821686957802283,
          -0.891604540812876101, -0.889489266395734002, -0.889489266395734002, -0.891604540812876101,
          -0.889821686957802283, -0.880029836830879542, -0.88270771158442396, -0.897476024922029492,
          -0.920665554037870759, -0.880560005368091314, -0.839125061489862656, -0.800413820399077491;
      }
      break;
      case 30:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << -0.800390085020429187, -0.839050608964517775, -0.880547772685447039, -0.920819039192321442,
          -0.897634652952583934, -0.882511895121084633, -0.879442500798295668, -0.889801872319090914,
          -0.892712487130870924, -0.890760710534785227, -0.887785410621690585, -0.888000317543117035,
          -0.889093931815407745, -0.889479357496824941, -0.889087453456860799, -0.889087453456860799,
          -0.88947935749682483, -0.889093931815407634, -0.888000317543117035, -0.887785410621690585,
          -0.890760710534785227, -0.892712487130871035, -0.889801872319090914, -0.879442500798295668,
          -0.882511895121084633, -0.897634652952583822, -0.920819039192321553, -0.880547772685447039,
          -0.839050608964517775, -0.800390085020429187;
      }
      break;
    }
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( -1 );
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo209 : public NonlinearSystem
{
public:
  RooseKullaLombMeressoo209( integer neq ) : NonlinearSystem( "Roose Kulla Lomb Meressoo N.209", RKM_BIBTEX, neq )
  {
    check_min_equations( n, 1 );
  }

  void checkx( Vector const & x ) const
  {
    for ( integer i = 0; i < n; ++i )
      UTILS_ASSERT( x( i ) > 0, "RooseKullaLombMeressoo209, found x[{}] = {} <= 0", i, x( i ) );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = power2( x( 0 ) ) - 1;
    for ( integer k = 1; k < n; ++k ) f( k ) = power2( x( k - 1 ) ) - 1 + log( x( k ) );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    J.insert( 0, 0 ) = 2 * x( 0 );
    for ( integer k = 1; k < n; ++k )
    {
      J.insert( k, k - 1 ) = 2 * x( k - 1 );
      J.insert( k, k )     = 1 / x( k );
    }
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 1 );
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
    for ( integer i = 1; i < n; ++i )
    {
      UTILS_ASSERT( x( i ) > 0, "Bad range" );
      UTILS_ASSERT( std::abs( x( i ) ) < 1000, "Bad range" );
    }
  }

  virtual void bounding_box( Vector & L, Vector & U ) const override
  {
    U.fill( 1000 );
    L.fill( -1000 );
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

// DA RICONTROLLARE!
class RooseKullaLombMeressoo210 : public NonlinearSystem
{
public:
  RooseKullaLombMeressoo210( integer neq ) : NonlinearSystem( "Roose Kulla Lomb Meressoo N.210", RKM_BIBTEX, neq )
  {
    check_min_equations( n, 1 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type acc = 0;
    for ( integer j = 0; j < n; ++j ) acc += x( j );
    for ( integer j = 0; j < n; ++j ) f( j ) = exp( cos( ( j + 1 ) * acc ) );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    real_type acc = 0;
    for ( integer j = 0; j < n; ++j ) acc += x( j );

    for ( integer i = 0; i < n; ++i )
    {
      real_type cc = cos( ( i + 1 ) * acc );
      real_type ss = sin( ( i + 1 ) * acc );
      for ( integer j = 0; j < n; ++j ) J.insert( i, j ) = -exp( cc ) * ss * ( i + 1 );
    }

    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.setZero();
  }

  string note() const { return "This problem has no solution"; }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo211 : public NonlinearSystem
{
public:
  RooseKullaLombMeressoo211( integer neq ) : NonlinearSystem( "Roose Kulla Lomb Meressoo N.211", RKM_BIBTEX, neq )
  {
    check_min_equations( n, 1 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type acc = 0;
    for ( integer j = 0; j < n; ++j ) acc += power3( x( j ) );
    for ( integer j = 0; j < n; ++j ) f( j ) = ( acc + j + 1 ) / ( 2 * n );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type tmp = 1.5 / n;
    for ( integer i = 0; i < n; ++i )
      for ( integer j = 0; j < n; ++j ) J.insert( i, j ) = tmp * power2( x( j ) );
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.setZero();
  }

  string note() const { return "This problem has no solution"; }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo212 : public NonlinearSystem
{
public:
  RooseKullaLombMeressoo212( integer neq ) : NonlinearSystem( "Roose Kulla Lomb Meressoo N.212", RKM_BIBTEX, neq )
  {
    check_min_equations( n, 1 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = x( 0 );
    for ( integer k = 1; k < n; ++k ) f( k ) = cos( x( k - 1 ) ) + x( k ) - 1;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    J.insert( 0, 0 ) = 1;
    for ( integer k = 1; k < n; ++k )
    {
      J.insert( k, k - 1 ) = -sin( x( k - 1 ) );
      J.insert( k, k )     = 1;
    }
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
    x0.fill( 0.5 );
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo213 : public NonlinearSystem
{
  real_type h2;

public:
  RooseKullaLombMeressoo213( integer neq ) : NonlinearSystem( "Roose Kulla Lomb Meressoo N.213", RKM_BIBTEX, neq )
  {
    check_min_equations( n, 1 );
    h2 = 1.0 / power2( n + 1 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    for ( integer k = 0; k < n; ++k )
    {
      real_type xp = k < n - 1 ? x( k + 1 ) : 1;
      real_type xm = k > 0 ? x( k - 1 ) : 0;
      real_type xc = x( k );
      f( k )       = 2 * xc - xm - xp + h2 * ( xc + sin( xc ) );
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer k = 0; k < n; ++k )
    {
      if ( k > 0 ) J.insert( k, k - 1 ) = -1;
      J.insert( k, k ) = 2 + h2 * ( 1 + cos( x( k ) ) );
      if ( k < n - 1 ) J.insert( k, k + 1 ) = -1;
    }
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    switch ( n )
    {
      case 2:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << 0.25493102645914878, 0.566207574177154838;
      }
      break;
      case 5:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << 0.123717884647979601, 0.254300224963985633, 0.398934465982092534, 0.565440128266567754,
          0.762535446654485805;
      }
      break;
      case 10:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << 0.0670031343864649215, 0.135113344358246179, 0.20545343835958349, 0.279177536667768345,
          0.35748628505027491, 0.441641370870779482, 0.532978814497350184, 0.632920225489453259, 0.74298082731276216,
          0.864772561262492578;
      }
      break;
      case 20:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << 0.0350196898973525322, 0.0701981830199509271, 0.105694904547656857, 0.141670522051372094,
          0.178287562727141746, 0.21571102540715037, 0.254108984809966743, 0.293653184781339449, 0.334519616353437577,
          0.376889075292977049, 0.42094769239294999, 0.466887428063573084, 0.514906520772035803, 0.565209876549352397,
          0.618009384118297755, 0.673524137209132956, 0.731980542357675512, 0.793612287002704475, 0.858660139154802993,
          0.927371546514711986;
      }
      break;
      case 30:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << 0.0237123159942052714, 0.0474739789290085645, 0.0713344245144125061, 0.095343265901786059,
          0.119550382155985915, 0.144006006415265148, 0.168760813611468746, 0.193866007602538554, 0.21937340754326351,
          0.245335533288195989, 0.271805689582330812, 0.298838048750088958, 0.32648773154088373, 0.354810885729579639,
          0.383864762001958348, 0.413707786578368653, 0.444399629942578145, 0.47600127094703848, 0.50857505546002324,
          0.542184748604255695, 0.576895579510812273, 0.612774277376699539, 0.649889097470416455, 0.688309835578470475,
          0.728107829229357573, 0.769355943873014314, 0.812128542037363621, 0.856501433334843409, 0.902551803057848212,
          0.950358116991871227;
      }
      break;
    }
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 1 );
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo214 : public NonlinearSystem
{
public:
  RooseKullaLombMeressoo214( integer neq ) : NonlinearSystem( "Roose Kulla Lomb Meressoo N.214", RKM_BIBTEX, neq )
  {
    check_min_equations( n, 1 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    integer zero = 0;
    for ( integer k = 0; k < n; ++k )
    {
      real_type ff   = 0;
      integer   imin = max( zero, k - 5 );
      integer   imax = min( n - 1, k + 1 );
      for ( integer i = imin; i <= imax; ++i )
        if ( i != k ) ff += x( i ) + power2( x( i ) );
      f( k ) = x( k ) * ( 2 + 5 * power2( x( k ) ) ) + 1 - ff;
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    integer zero = 0;
    for ( integer k = 0; k < n; ++k )
    {
      integer imin = max( zero, k - 5 );
      integer imax = min( n - 1, k + 1 );
      for ( integer i = imin; i <= imax; ++i )
      {
        if ( i == k )
          J.insert( k, i ) = 2 + 15 * power2( x( k ) );
        else
          J.insert( k, i ) = -( 1 + 2 * x( i ) );
      }
    }
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    switch ( n )
    {
      case 2:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << -0.427304623558166286, -0.427304623558166286;
      }
      break;
      case 5:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << -0.428302864642700787, -0.476596531501095377, -0.519637722100754651, -0.558861956527025194,
          -0.558861956527025194;
      }
      break;
      case 10:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << -0.428302863587250282, -0.476596424356290238, -0.519652463646861684, -0.558099324832180832,
          -0.592506156829457287, -0.624503682199467947, -0.623239471440591108, -0.621393841796573532,
          -0.620453596659087392, -0.586469270720434976;
      }
      break;
      case 20:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << -0.428302863587250282, -0.476596424356293569, -0.519652463646401386, -0.55809932485615199,
          -0.592506155965082826, -0.624503707410516529, -0.623238669132451184, -0.621419676713647839,
          -0.61961584283347626, -0.618226017919857429, -0.617518024841495761, -0.617731830318644759,
          -0.617900316253351289, -0.618007798540867848, -0.618057061755049264, -0.618062699716298014,
          -0.618047199350808651, -0.618011195738616514, -0.618872079495047522, -0.586276945400115101;
      }
      break;
      case 30:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << -0.428302863587250282, -0.476596424356293569, -0.519652463646401386, -0.55809932485615199,
          -0.592506155965082826, -0.624503707410516529, -0.623238669132451184, -0.621419676713647839,
          -0.61961584283347626, -0.618226017919857429, -0.617518024841495206, -0.617731830318665742,
          -0.617900316252663617, -0.61800779856335919, -0.618057061019478993, -0.618062723774471579,
          -0.61804641236762925, -0.618036943255954929, -0.618032796823900332, -0.618032010907616169,
          -0.618032748437421176, -0.61803365220978157, -0.618034039196207474, -0.618034129052205672,
          -0.618034091025163379, -0.618034003909173957, -0.618034776213912562, -0.618008230615912701,
          -0.618873272626757731, -0.58627911806458255;
      }
      break;
    }
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( -1 );
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo215 : public NonlinearSystem
{
  mutable Vector grad_dS;
  mutable Vector grad_S;

  using Matrix = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;

public:
  RooseKullaLombMeressoo215( integer neq ) : NonlinearSystem( "Roose Kulla Lomb Meressoo N.215", RKM_BIBTEX, neq )
  {
    check_min_equations( n, 1 );
    grad_dS.resize( n );
    grad_S.resize( n );
  }

  real_type g( integer k ) const { return ( k + 1.0 ) / 29.0; }

  real_type S( Vector const & x, integer k ) const
  {
    real_type sum = 0;
    real_type gkj = 1;
    real_type gk  = g( k );
    for ( integer j = 0; j < n; ++j, gkj *= gk ) sum += gkj * x( j );
    return sum;
  }

  void S_1( Vector const &, integer k, Vector & grad ) const
  {
    real_type gkj = 1;
    real_type gk  = g( k );
    for ( integer j = 0; j < n; ++j, gkj *= gk ) grad[j] = gkj;
  }

  real_type dS( Vector const & x, integer k ) const
  {
    real_type sum = 0;
    real_type gkj = 1;
    real_type gk  = g( k );
    for ( integer j = 1; j < n; ++j, gkj *= gk ) sum += j * gkj * x( j );
    return sum;
  }

  void dS_1( Vector const &, integer k, Vector & grad ) const
  {
    real_type gkj = 1;
    real_type gk  = g( k );
    grad[0]       = 0;
    for ( integer j = 1; j < n; ++j, gkj *= gk ) grad[j] = j * gkj;
  }

  real_type powergk( integer k, integer i ) const
  {
    real_type gk = g( k );
    if ( i == 0 ) return 1 / gk;
    real_type res = 1;
    for ( integer j = 1; j < i; ++j ) res *= gk;
    return res;
  }

  virtual void evaluate( Vector const & x, Vector & ff ) const override
  {
    for ( integer i = 0; i < n; ++i )
    {
      real_type f = 0;
      for ( integer k = 0; k < 29; ++k )
      {
        real_type Sk    = S( x, k );
        real_type dSk   = dS( x, k );
        real_type gk    = g( k );
        real_type powgk = powergk( k, i );
        f += powgk * ( i - 2 * gk * Sk ) * ( dSk - Sk * Sk - 1 );
      }
      if ( i == 0 ) f += x( 0 ) * ( 1 - 2 * ( x( 1 ) - x( 0 ) * x( 0 ) - 1 ) );
      if ( i == 1 ) f += x( 1 ) - x( 0 ) * x( 0 ) - 1;
      ff( i ) = f;
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    Matrix J_full( n, n );
    J_full.setZero();
    for ( integer i = 0; i < n; ++i )
    {
      for ( integer k = 0; k < 29; ++k )
      {
        real_type gk    = g( k );
        real_type powgk = powergk( k, i );
        real_type Sk    = S( x, k );
        real_type dSk   = dS( x, k );
        dS_1( x, k, grad_dS );
        S_1( x, k, grad_S );
        real_type A = i - 2 * gk * Sk;
        real_type B = dSk - Sk * Sk - 1;
        for ( integer j = 0; j < n; ++j )
        {
          real_type A_1 = -2 * gk * grad_S[j];
          real_type B_1 = grad_dS[j] - 2 * Sk * grad_S[j];
          J_full( i, j ) += powgk * ( A_1 * B + A * B_1 );
        }
      }
    }
    J_full( 0, 0 ) += 6 * x( 0 ) * x( 0 ) - 2 * x( 1 ) + 3;
    J_full( 0, 1 ) -= 2 * x( 0 );
    J_full( 1, 0 ) -= 2 * x( 0 );
    J_full( 1, 1 ) += 1;
    J.resize( n, n );
    J = J_full.sparseView();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    switch ( n )
    {
      case 2:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << -0.501367007521962726, 1.0736498384595432;
      }
      break;
      case 5:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << -0.0714051701051878485, 0.970024704720390041, 0.2661688465197945, -0.539893187056864843,
          0.710711221847347585;
      }
      break;
      case 6:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << -0.0157250864014588515, 1.01243486936910987, -0.232991625956739556, 1.26043008779961352,
          -1.51372892272228565, 0.992996432431136333;
      }
      break;
      case 10:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << -1.22248986826704962e-06, 1.00004450564633962, -0.00508490059336646257, 0.417429340147828398,
          -0.588443484856691734, 2.30197280176861696, -4.56245409953258463, 5.59197371897036266, -3.64040714755196371,
          1.04237104657441626;
      }
      break;
    }
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.setZero();
    // get_exact_solution( x, 0 );
  }

  string note() const { return "Watson Function"; }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo216 : public NonlinearSystem
{
  mutable std::vector<std::vector<real_type>> m_y, m_y_D;

public:
  RooseKullaLombMeressoo216( integer neq ) : NonlinearSystem( "Roose Kulla Lomb Meressoo N.216", RKM_BIBTEX, neq )
  {
    UTILS_ASSERT( n > 0 && n != 8 && n < 10, "RooseKullaLombMeressoo216, neq={} must be [1..7] or 9", n );
    m_y.resize( n );
    m_y_D.resize( n );
    for ( integer i = 0; i < n; ++i )
    {
      m_y[i].resize( n );
      m_y_D[i].resize( n );
    }
  }

  void evalY( real_type s, std::vector<real_type> & y ) const
  {
    y[0] = 2 * s - 1;
    y[1] = 8 * ( s - 1 ) * s + 1;
    for ( integer k = 2; k < n; ++k ) y[k] = 2 * y[0] * y[k - 1] - y[k - 2];
  }

  void evalY_D( real_type s, std::vector<real_type> & y, std::vector<real_type> & y_D ) const
  {
    evalY( s, y );
    y_D[0] = 2;
    y_D[1] = 16 * s - 8;
    for ( integer k = 2; k < n; ++k ) y_D[k] = 2 * ( y[0] * y_D[k - 1] + y_D[0] * y[k - 1] ) - y_D[k - 2];
  }

  void evalY( Vector const & x ) const
  {
    for ( integer k = 0; k < n; ++k ) evalY( x( k ), m_y[k] );
  }

  void evalY_D( Vector const & x ) const
  {
    for ( integer k = 0; k < n; ++k ) evalY_D( x( k ), m_y[k], m_y_D[k] );
  }

  real_type Y( integer xj, integer k ) const { return m_y[xj][k]; }

  real_type Y_D( integer xj, integer k ) const { return m_y_D[xj][k]; }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    evalY( x );
    for ( integer k = 0; k < n; ++k )
    {
      f( k ) = 0;
      for ( integer j = 0; j < n; ++j ) f( k ) += Y( j, k );
      f( k ) /= n;
      if ( ( k % 2 ) == 1 ) f( k ) += 1.0 / ( power2( k + 1 ) - 1.0 );
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    evalY_D( x );
    for ( integer k = 0; k < n; ++k )
      for ( integer j = 0; j < n; ++j ) J.insert( k, j ) = Y_D( j, k ) / n;
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    switch ( n )
    {
      case 2:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << 0.211324865405187134, 0.788675134594812866;
      }
      break;
      case 3:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << 0.146446609406726241, 0.5, 0.853553390593273731;
      }
      break;
      case 4:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << 0.10267276385411693, 0.406203762957460079, 0.593796237042539921, 0.897327236145883056;
      }
      break;
      case 5:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << 0.0837512564995090553, 0.312729295223209469, 0.5, 0.687270704776790531, 0.916248743500490903;
      }
      break;
      case 6:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << 0.0668765909460896923, 0.288740673119444125, 0.366682299241647691, 0.633317700758352253,
          0.71125932688055582, 0.933123409053910335;
      }
      break;
      case 7:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << 0.0580691496209754937, 0.23517161235742165, 0.338044094740046208, 0.5, 0.661955905259953847,
          0.764828387642578433, 0.941930850379024465;
      }
      break;
      case 9:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << 0.0442053461357827734, 0.199490672309881045, 0.235619108471059907, 0.416046907892598183,
          0.499999999999999833, 0.583953092107402094, 0.764380891528939954, 0.800509327690119066, 0.955794653864217247;
      }
      break;
    }
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    real_type bf = 1.0 / ( n + 1 );
    for ( integer k{ 0 }; k < n; ++k ) x0( k ) = bf * ( k + 1 );
  }

  virtual void check_if_admissible( Vector const & x ) const override
  {
    for ( integer i = 0; i < n; ++i ) UTILS_ASSERT( x( i ) > 0, "Bad range" );
  }

  string note() const { return "Each permutation of x is a solution"; }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo217 : public NonlinearSystem
{
public:
  RooseKullaLombMeressoo217( integer neq ) : NonlinearSystem( "Roose Kulla Lomb Meressoo N.217", RKM_BIBTEX, neq )
  {
    check_min_equations( n, 2 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    for ( integer k = 0; k < n; ++k )
    {
      real_type xc = x( k );
      real_type xp = k < n - 1 ? x( k + 1 ) : 20;
      real_type xm = k > 0 ? x( k - 1 ) : 0;
      f( k )       = 3 * xc * ( xm + xp - 2 * xc ) + power2( xp - xm ) / 4;
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer k = 0; k < n; ++k )
    {
      real_type xc     = x( k );
      real_type xp     = k < n - 1 ? x( k + 1 ) : 20;
      real_type xm     = k > 0 ? x( k - 1 ) : 0;
      J.insert( k, k ) = 3 * ( xp + xm - 4 * xc );
      if ( k > 0 ) J.insert( k, k - 1 ) = 3 * xc - ( xp - xm ) / 2;
      if ( k < n - 1 ) J.insert( k, k + 1 ) = 3 * xc + ( xp - xm ) / 2;
    }
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    switch ( n )
    {
      case 2:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << 8.33828484333201025, 14.5583676083251294;
      }
      break;
      case 5:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << 4.88928869829531543, 8.53653521682389993, 11.7273247200356803, 14.6523320192791218, 17.394667732103354;
      }
      break;
      case 10:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << 3.08315248959636978, 5.38308155447113368, 7.39517190291697979, 9.23966178544206862, 10.9689601971417012,
          12.6118651601462197, 14.1863707080987851, 15.7046865038081975, 17.1755885168751945, 18.6056591191924703;
      }
      break;
      case 20:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << 1.89123927553494386, 3.30204078247077248, 4.53627888965039894, 5.66770904788268393, 6.72847950486216018,
          7.73625528062737722, 8.70207411115346829, 9.63342553642360322, 10.53569283165554, 11.4129136953710457,
          12.2682175592459135, 13.1040939164598367, 13.9225656143308374, 14.7253054956124636, 15.5137176336911349,
          16.2889955537418949, 17.0521649909861708, 17.8041159610068647, 18.5456272591407938, 19.2773854806811933;
      }
      break;
      case 30:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << 1.41029144531941619, 2.46232189012527636, 3.38268953823807861, 4.22639360771921346, 5.01740695028136763,
          5.76890231829069311, 6.48911051829969399, 7.18361647777063173, 7.85643448885871631, 8.51057545121989989,
          9.14837297265207816, 9.77168346073938032, 10.3820153466489788, 10.9806160641959405, 11.5685326266858191,
          12.1466550422555599, 12.7157481985778418, 13.2764757775986784, 13.8294185246607579, 14.3750884318112888,
          14.9139399076884889, 15.4463786872421789, 15.9727690205710928, 16.4934395336305428, 17.0086880512958594,
          17.5187856006541267, 18.0239797600299596, 18.5244974809442624, 19.0205474818159175, 19.5123222909239296;
      }
      break;
    }
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 10 );
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo218 : public NonlinearSystem
{
  real_type h, h2;

public:
  RooseKullaLombMeressoo218( integer neq ) : NonlinearSystem( "Roose Kulla Lomb Meressoo N.218", RKM_BIBTEX, neq )
  {
    check_min_equations( n, 2 );
    h  = 1.0 / ( n + 1 );
    h2 = h * h;
  }

  real_type t( integer j ) const { return ( j + 1 ) * h; }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    for ( integer k = 0; k < n; ++k )
    {
      real_type xc = x( k );
      real_type xp = k < n - 1 ? x( k + 1 ) : 0;
      real_type xm = k > 0 ? x( k - 1 ) : 0;
      f( k )       = 2 * xc - ( xm + xp ) + ( 0.5 * h2 ) * power3( xc + t( k ) + 1 );
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer k = 1; k < n; ++k ) J.insert( k, k - 1 ) = -1;
    for ( integer k = 0; k < n - 1; ++k ) J.insert( k, k + 1 ) = -1;
    for ( integer k = 0; k < n; ++k ) J.insert( k, k ) = 2 + 1.5 * h2 * power2( x( k ) + t( k ) + 1 );
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    switch ( n )
    {
      case 2:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << -0.128246763033731614, -0.159267567244640834;
      }
      break;
      case 5:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << -0.0750221292923204525, -0.131976210352190593, -0.164848771909337249, -0.164664680215800663,
          -0.117417651684193575;
      }
      break;
      case 10:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << -0.0431649825187648689, -0.0815771565353868716, -0.114485714380529277, -0.14097357686259665,
          -0.159908696181983112, -0.169877202312774922, -0.169089983781208347, -0.155249535221831853,
          -0.125355891678934989, -0.0754165336858920871;
      }
      break;
      case 20:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << -0.0232104399966753146, -0.0452020277196446968, -0.0658809801847738824, -0.0851436508063508207,
          -0.102875197721159828, -0.118948030435962207, -0.133219990530385718, -0.145532211750480533,
          -0.155706591602235406, -0.16354278961733093, -0.16881464562153628, -0.171265882959322213,
          -0.170604924470161234, -0.166498599947109888, -0.158564458442799261, -0.146361310877112738,
          -0.129377508964773374, -0.107016302444755948, -0.0785773886624082651, -0.0432334478445528594;
      }
      break;
      case 30:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << -0.0158588747608703062, -0.0311714390223493919, -0.0459099102817520852, -0.060044590713602998,
          -0.073543699225747064, -0.0863731855306687224, -0.0984965239444887258, -0.10987448428747057,
          -0.120464876863777437, -0.130222268033639288, -0.139097662344602108, -0.147038146543775367,
          -0.153986490029936474, -0.159880695398378264, -0.164653491652127226, -0.168231761362994642,
          -0.170535891518048321, -0.17147903592302513, -0.17096627478051718, -0.168893654324848153,
          -0.165147086059982073, -0.159601081061931882, -0.15211728978118419, -0.142542811566464672,
          -0.130708230408252885, -0.116425323750638757, -0.0994843790941261352, -0.0796510377833258565,
          -0.0566625658741519017, -0.0302234270054019573;
      }
      break;
    }
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    for ( integer k{ 0 }; k < n; ++k ) x0( k ) = t( k ) * ( t( k ) - 1 );
  }

  string note() const { return "For any n all solutions are in the interval [-0.5,0]"; }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo219 : public NonlinearSystem
{
  real_type h;

  using Matrix = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;

public:
  RooseKullaLombMeressoo219( integer neq ) : NonlinearSystem( "Roose Kulla Lomb Meressoo N.219", RKM_BIBTEX, neq )
  {
    check_min_equations( n, 2 );
    h = 1.0 / ( n + 1 );
  }

  real_type t( integer j ) const { return ( j + 1 ) * h; }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    for ( integer i = 0; i < n; ++i )
    {
      real_type sum1 = 0;
      real_type sum2 = 0;
      for ( integer j = 0; j <= i; ++j ) sum1 += t( j ) * power3( x( j ) + t( j ) + 1 );
      for ( integer j = i + 1; j < n; ++j ) sum2 += ( 1 - t( j ) ) * power3( x( j ) + t( j ) + 1 );
      f( i ) = x( i ) + 0.5 * h * ( ( 1 - t( i ) ) * sum1 + t( i ) * sum2 );
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    Matrix J_full( n, n );
    J_full.setZero();

    for ( integer i = 0; i < n; ++i )
    {
      real_type bf = 1.5 * h * ( 1 - t( i ) );
      for ( integer j = 0; j <= i; ++j ) J_full( i, j ) = bf * t( j ) * power2( x( j ) + t( j ) + 1 );

      J_full( i, i ) += 1;  // somma su j(i,i)

      bf = 1.5 * h * t( i );
      for ( integer j = i + 1; j < n; ++j ) J_full( i, j ) = bf * ( 1 - t( j ) ) * power2( x( j ) + t( j ) + 1 );
    }

    J.resize( n, n );
    J = J_full.sparseView();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    switch ( n )
    {
      case 2:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << -0.128246763033731614, -0.159267567244640834;
      }
      break;
      case 5:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << -0.0750221292923204663, -0.13197621035219062, -0.164848771909337277, -0.164664680215800718,
          -0.117417651684193616;
      }
      break;
      case 10:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << -0.043164982518764862, -0.0815771565353868855, -0.114485714380529277, -0.14097357686259665,
          -0.159908696181983112, -0.169877202312774894, -0.169089983781208347, -0.155249535221831825,
          -0.125355891678934933, -0.0754165336858920177;
      }
      break;
      case 20:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << -0.0232104399966752729, -0.0452020277196446066, -0.0658809801847737575, -0.0851436508063506542,
          -0.102875197721159606, -0.118948030435961943, -0.133219990530385413, -0.145532211750480228,
          -0.1557065916022351, -0.163542789617330597, -0.168814645621535919, -0.171265882959321825,
          -0.170604924470160901, -0.166498599947109527, -0.158564458442798956, -0.146361310877112488,
          -0.129377508964773152, -0.107016302444755837, -0.0785773886624082513, -0.0432334478445528664;
      }
      break;
      case 30:
      {
        x_vec.resize( 1 );
        auto & x0{ x_vec[0] };
        x0.resize( n );
        x0 << -0.0158588747608703166, -0.0311714390223494162, -0.0459099102817521199, -0.0600445907136030396,
          -0.0735436992257471334, -0.0863731855306687918, -0.0984965239444888091, -0.109874484287470653,
          -0.120464876863777492, -0.130222268033639343, -0.139097662344602191, -0.147038146543775394,
          -0.153986490029936474, -0.159880695398378181, -0.164653491652127115, -0.168231761362994447,
          -0.170535891518048127, -0.171479035923024881, -0.170966274780516903, -0.168893654324847792,
          -0.165147086059981713, -0.159601081061931466, -0.152117289781183801, -0.142542811566464256,
          -0.130708230408252524, -0.116425323750638438, -0.0994843790941258299, -0.0796510377833256483,
          -0.0566625658741517907, -0.0302234270054018567;
      }
      break;
    }
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    for ( integer k{ 0 }; k < n; ++k ) x0( k ) = t( k ) * ( t( k ) - 1 );
  }
};
