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

class EsterificReaction : public NonlinearSystem
{
public:
  EsterificReaction()
    : NonlinearSystem(
        "EsterificReaction (example 4)",
        "@book{Luus:2000,\n"
        "  author    = {Luus, Rein},\n"
        "  title     = {Iterative Dynamic Programming},\n"
        "  year      = {2000},\n"
        "  isbn      = {1584881488},\n"
        "  edition   = {1st},\n"
        "  publisher = {CRC Press, Inc.},\n"
        "}\n",
        2 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type X   = x( 0 );
    real_type Y   = x( 1 );
    real_type t1  = 0.125E1 * X;
    real_type t2  = 0.225E1 * Y;
    real_type t12 = ( X + Y ) * ( 0.6875E1 * X + 0.7875E1 * Y ) / ( 0.5797E1 * X + 0.6797E1 * Y );
    real_type t16 = pow( t12 - Y - X, 2.0 );
    real_type t20 = pow( 0.11526E1 * t12 - t1 - t2, 2.0 );
    f( 0 )        = ( t1 + t2 - 0.1054E1 * t12 ) * t16 - 0.20564E4 * t20;
    f( 1 )        = 0.55E1 * ( -0.1433314531E3 * t12 + 0.2801659319E3 * Y + 0.7647126E2 + 0.1556548821E3 * X ) *
               ( -0.1043629412E3 * t12 + 0.2031305325E3 * Y + 0.1134099326E3 * X + 0.61177E2 ) -
             ( 0.1422774531E3 * t12 - 0.2789159319E3 * Y - 0.1544048821E3 * X - 0.7585949E2 ) *
               ( 0.11526E3 * t12 - 0.225E3 * Y - 0.125E3 * X - 0.61177E2 );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type X      = x( 0 );
    real_type Y      = x( 1 );
    real_type t3     = 0.6875E1 * X + 0.7875E1 * Y;
    real_type t6     = 0.5797E1 * X + 0.6797E1 * Y;
    real_type t7     = 1 / t6;
    real_type t8     = t3 * t7;
    real_type t9     = 0.1054E1 * t8;
    real_type t10    = X + Y;
    real_type t11    = t10 * t7;
    real_type t13    = t10 * t3;
    real_type t14    = t6 * t6;
    real_type t16    = t13 / t14;
    real_type t19    = t13 * t7;
    real_type t20    = t19 - Y - X;
    real_type t21    = t20 * t20;
    real_type t23    = 0.125E1 * X;
    real_type t24    = 0.225E1 * Y;
    real_type t27    = ( t23 + t24 - 0.1054E1 * t19 ) * t20;
    real_type t34    = 0.11526E1 * t19 - t23 - t24;
    real_type t35    = 0.11526E1 * t8;
    real_type t57    = 0.1433314531E3 * t8;
    real_type t64    = -0.1043629412E3 * t19 + 0.2031305325E3 * Y + 0.1134099326E3 * X + 0.61177E2;
    real_type t70    = -0.1433314531E3 * t19 + 0.2801659319E3 * Y + 0.7647126E2 + 0.1556548821E3 * X;
    real_type t71    = 0.1043629412E3 * t8;
    real_type t77    = 0.1422774531E3 * t8;
    real_type t84    = 0.11526E3 * t19 - 0.225E3 * Y - 0.125E3 * X - 0.61177E2;
    real_type t89    = 0.1422774531E3 * t19 - 0.2789159319E3 * Y - 0.1544048821E3 * X - 0.7585949E2;
    real_type t90    = 0.11526E3 * t8;
    J.insert( 0, 0 ) = ( 0.125E1 - t9 - 0.724625E1 * t11 + 0.6110038E1 * t16 ) * t21 +
                       2.0 * t27 * ( t8 + 0.6875E1 * t11 - 0.5797E1 * t16 - 1.0 ) -
                       0.41128E4 * t34 * ( t35 + 0.7924125E1 * t11 - 0.66816222E1 * t16 - 0.125E1 );
    J.insert( 0, 1 ) = ( 0.225E1 - t9 - 0.830025E1 * t11 + 0.7164038E1 * t16 ) * t21 +
                       2.0 * t27 * ( t8 + 0.7875E1 * t11 - 0.6797E1 * t16 - 1.0 ) -
                       0.41128E4 * t34 * ( t35 + 0.9076725E1 * t11 - 0.78342222E1 * t16 - 0.225E1 );
    J.insert( 1, 0 ) = 0.55E1 * ( -t57 - 0.9854037401E3 * t11 + 0.8308924336E3 * t16 + 0.1556548821E3 ) * t64 +
                       0.55E1 * t70 * ( -t71 - 0.7174952208E3 * t11 + 0.6049919701E3 * t16 + 0.1134099326E3 ) -
                       ( t77 + 0.9781574901E3 * t11 - 0.8247823956E3 * t16 - 0.1544048821E3 ) * t84 -
                       t89 * ( t90 + 0.7924125E3 * t11 - 0.66816222E3 * t16 - 0.125E3 );
    J.insert( 1, 1 ) = 0.55E1 * ( -t57 - 0.1128735193E4 * t11 + 0.9742238867E3 * t16 + 0.2801659319E3 ) * t64 +
                       0.55E1 * t70 * ( -t71 - 0.821858162E3 * t11 + 0.7093549113E3 * t16 + 0.2031305325E3 ) -
                       ( t77 + 0.1120434943E4 * t11 - 0.9670598487E3 * t16 - 0.2789159319E3 ) * t84 -
                       t89 * ( t90 + 0.9076725E3 * t11 - 0.78342222E3 * t16 - 0.225E3 );
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 3 );
    auto & x0{ x_vec[0] };
    auto & x1{ x_vec[1] };
    auto & x2{ x_vec[2] };
    x0.resize( n );
    x1.resize( n );
    x2.resize( n );
    x0 << 61.15373203411778822477588185225287963418, 6.976291974526651896957787060226279608637;
    x1 << 67.09266776655127299522316830183801702230, 7.614805557040459414336978005675721174002;
    x2 << -33.28012846942180982310942345376518595181, 28.52767665376609064013743556335143862501;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 20 );
  }
};
