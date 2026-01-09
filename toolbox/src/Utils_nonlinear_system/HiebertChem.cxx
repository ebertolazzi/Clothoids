/*\
 |
 |  Author:
 |    Enrico Bertolazzi
 |    University of Trento
 |    Department of Industrial Engineering
 |    Via Sommarive 9, I-38123, Povo, Trento, Italy
 |    email: enrico.bertolazzi@unitn.it
\*/

#define HIEBERT_BIBTEX                                         \
  "@article{Hiebert:1982,\n"                                   \
  "  author  = {Hiebert, K. L.},\n"                            \
  "  title   = {An Evaluation of Mathematical Software That\n" \
  "             Solves Systems of Nonlinear Equations},\n"     \
  "  journal = {ACM Trans. Math. Softw.},\n"                   \
  "  year    = {1982},\n"                                      \
  "  volume  = {8},\n"                                         \
  "  number  = {1},\n"                                         \
  "  pages   = {5--20},\n"                                     \
  "  doi     = {10.1145/355984.355986},\n"                     \
  "}\n"

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class HiebertChem2x2 : public NonlinearSystem
{
public:
  HiebertChem2x2() : NonlinearSystem( "Hiebert Chem 2x2", HIEBERT_BIBTEX, 2 ) {}

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = x( 1 ) - 10;
    f( 1 ) = x( 0 ) * x( 1 ) - 5E4;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    J.insert( 0, 0 ) = 0;
    J.insert( 0, 1 ) = 1;
    J.insert( 1, 0 ) = x( 1 );
    J.insert( 1, 1 ) = x( 0 );
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 5000, 10;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 1, 1;
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class HiebertChem6x6 : public NonlinearSystem
{
public:
  HiebertChem6x6() : NonlinearSystem( "Hiebert Chem 6x6", HIEBERT_BIBTEX, 6 ) {}

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = x( 0 ) + x( 1 ) + x( 3 ) - 0.001;
    f( 1 ) = x( 4 ) + x( 5 ) - 55;
    f( 2 ) = x( 0 ) + x( 1 ) + x( 2 ) + 2 * x( 4 ) + x( 5 ) - 110.001;
    f( 3 ) = x( 0 ) - 0.1 * x( 1 );
    f( 4 ) = x( 0 ) - 10000 * x( 2 ) * x( 3 );
    f( 5 ) = x( 3 ) - 55E14 * x( 2 ) * x( 5 );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    // Costanti sicure in floating point
    constexpr real_type C1 = 10000.0;
    constexpr real_type C2 = 5.5e15;

    // -----------------------------
    // f0 = x0 + x1 + x3 - 0.001
    // -----------------------------
    J.insert( 0, 0 ) = 1.0;
    J.insert( 0, 1 ) = 1.0;
    J.insert( 0, 3 ) = 1.0;

    // -----------------------------
    // f1 = x4 + x5 - 55
    // -----------------------------
    J.insert( 1, 4 ) = 1.0;
    J.insert( 1, 5 ) = 1.0;

    // -----------------------------
    // f2 = x0 + x1 + x2 + 2*x4 + x5 - 110.001
    // -----------------------------
    J.insert( 2, 0 ) = 1.0;
    J.insert( 2, 1 ) = 1.0;
    J.insert( 2, 2 ) = 1.0;
    J.insert( 2, 4 ) = 2.0;
    J.insert( 2, 5 ) = 1.0;

    // -----------------------------
    // f3 = x0 - 0.1*x1
    // -----------------------------
    J.insert( 3, 0 ) = 1.0;
    J.insert( 3, 1 ) = -0.1;

    // -----------------------------
    // f4 = x0 - 10000 * x2 * x3
    // -----------------------------
    real_type x2_safe = std::abs( x( 2 ) );
    real_type x3_safe = std::abs( x( 3 ) );
    J.insert( 4, 0 )  = 1.0;
    J.insert( 4, 2 )  = -C1 * x3_safe;
    J.insert( 4, 3 )  = -C1 * x2_safe;

    // -----------------------------
    // f5 = x3 - 5.5e15 * x2 * x5
    // -----------------------------
    real_type x2_f5  = std::abs( x( 2 ) );
    real_type x5_f5  = std::abs( x( 5 ) );
    J.insert( 5, 2 ) = -C2 * x5_f5;
    J.insert( 5, 3 ) = 1.0;
    J.insert( 5, 5 ) = -C2 * x2_f5;

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

    x0 << 0.10000000000000865800865802832637712049655126051610e-3,
      0.10000000000000865800865802832637712049655126051610e-2, -0.99999999999913419913419799193454109698763804853664e-4,
      -0.10000000000009523809523831159014832546206386567714e-3, 54.999999999999999818181818181487603305784236699939,
      0.18181818181851239669421576330006082348099514876998e-15;

    x1 << -0.18181818182178512396701422689707131359333381750224e-14,
      -0.18181818182178512396701422689707131359333381750224e-13,
      -0.18181818181814876033057918557475582208319103887114e-15,
      0.10000000000200000000003963636363715649586778444953e-2, 55.001000000000020181818182214512396702144144252600,
      -0.10000000000201818181822145123967021441442526003173e-2;

    x2 << 0.82644628099181424635970070847176598631923570438011e-4,
      0.82644628099181424635970070847176598631923570438011e-3, 0.90909090909186147186147038862875597975247651584838e-4,
      0.90909090909004329004329220681057415048840725181878e-4, 54.999999999999999818181818182181818181817073593074,
      0.18181818181781818181818292640692640295982173534745e-15;
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

class HiebertChem10x10 : public NonlinearSystem
{
  real_type const R;

public:
  HiebertChem10x10() : NonlinearSystem( "Hiebert Chem 10x10", HIEBERT_BIBTEX, 10 ), R( 3 ) {}

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type TOT = x( 0 ) + x( 1 ) + x( 2 ) + x( 3 ) + x( 4 ) + x( 5 ) + x( 6 ) + x( 7 ) + x( 8 ) + x( 9 );
    f( 0 )        = x( 0 ) + x( 3 ) - 3.0;
    f( 1 )        = 2.0 * x( 0 ) + x( 1 ) + x( 3 ) + x( 6 ) + x( 7 ) + x( 8 ) + 2.0 * x( 9 ) - R;
    f( 2 )        = 2.0 * x( 1 ) + 2.0 * x( 4 ) + x( 5 ) + x( 6 ) - 8.0;
    f( 3 )        = 2.0 * x( 2 ) + x( 4 ) - 4.0 * R;
    f( 4 )        = x( 0 ) * x( 4 ) - ( 193.0 / 1000.0 ) * x( 1 ) * x( 3 );
    f( 5 )        = x( 5 ) * sqrt( x( 1 ) ) - ( 2597.0 / 1000000.0 ) * sqrt( x( 1 ) * x( 3 ) * TOT );
    f( 6 )        = x( 6 ) * sqrt( x( 3 ) ) - ( 431.0 / 125000.0 ) * sqrt( x( 0 ) * x( 3 ) * TOT );
    f( 7 )        = x( 7 ) * x( 3 ) - ( 1799.0 / 100000000.0 ) * x( 2 ) * TOT;
    f( 8 )        = x( 8 ) * x( 3 ) - ( 431.0 / 2000000.0 ) * x( 0 ) * sqrt( x( 2 ) * TOT );
    f( 9 )        = ( x( 9 ) - ( 1923.0 / 50000000.0 ) * TOT ) * x( 3 ) * x( 3 );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    real_type t3  = x( 1 ) * x( 3 );
    real_type t4  = x( 0 ) + x( 1 ) + x( 2 ) + x( 3 ) + x( 4 ) + x( 5 ) + x( 6 ) + x( 7 ) + x( 8 ) + x( 9 );
    real_type t6  = sqrt( t3 * t4 );
    real_type t7  = 1 / t6;
    real_type t10 = 2597.0 / 2000000.0 * t7 * x( 1 ) * x( 3 );
    real_type t11 = sqrt( x( 1 ) );
    real_type t15 = x( 3 ) * t4;
    real_type t25 = x( 0 ) * x( 3 );
    real_type t27 = sqrt( t25 * t4 );
    real_type t28 = 1 / t27;
    real_type t34 = 431.0 / 250000.0 * t28 * x( 0 ) * x( 3 );
    real_type t35 = sqrt( x( 3 ) );
    real_type t45 = 1799.0 / 100000000.0 * x( 2 );
    real_type t60 = sqrt( x( 2 ) * t4 );
    real_type t63 = x( 0 ) / t60;
    real_type t65 = 431.0 / 4000000.0 * t63 * x( 2 );
    real_type t73 = x( 3 ) * x( 3 );
    real_type t74 = 1923.0 / 50000000.0 * t73;

    J.insert( 0, 0 ) = 1.0;
    J.insert( 0, 3 ) = 1.0;
    J.insert( 1, 0 ) = 2.0;
    J.insert( 1, 1 ) = 1.0;
    J.insert( 1, 3 ) = 1.0;
    J.insert( 1, 6 ) = 1.0;
    J.insert( 1, 7 ) = 1.0;
    J.insert( 1, 8 ) = 1.0;
    J.insert( 1, 9 ) = 2.0;
    J.insert( 2, 1 ) = 2.0;
    J.insert( 2, 4 ) = 2.0;
    J.insert( 2, 5 ) = 1.0;
    J.insert( 2, 6 ) = 1.0;
    J.insert( 3, 2 ) = 2.0;
    J.insert( 3, 4 ) = 1.0;
    J.insert( 4, 0 ) = x( 4 );
    J.insert( 4, 1 ) = -193.0 / 1000.0 * x( 3 );
    J.insert( 4, 3 ) = -193.0 / 1000.0 * x( 1 );
    J.insert( 4, 4 ) = x( 0 );
    J.insert( 5, 0 ) = -t10;
    J.insert( 5, 1 ) = x( 5 ) / t11 / 2.0 - 2597.0 / 2000000.0 * t7 * ( t15 + t3 );
    J.insert( 5, 2 ) = -t10;
    J.insert( 5, 3 ) = -2597.0 / 2000000.0 * t7 * ( x( 1 ) * t4 + t3 );
    J.insert( 5, 4 ) = -t10;
    J.insert( 5, 5 ) = t11 - t10;
    J.insert( 5, 6 ) = -t10;
    J.insert( 5, 7 ) = -t10;
    J.insert( 5, 8 ) = -t10;
    J.insert( 5, 9 ) = -t10;
    J.insert( 6, 0 ) = -431.0 / 250000.0 * t28 * ( t15 + t25 );
    J.insert( 6, 1 ) = -t34;
    J.insert( 6, 2 ) = -t34;
    J.insert( 6, 3 ) = x( 6 ) / t35 / 2.0 - 431.0 / 250000.0 * t28 * ( x( 0 ) * t4 + t25 );
    J.insert( 6, 4 ) = -t34;
    J.insert( 6, 5 ) = -t34;
    J.insert( 6, 6 ) = t35 - t34;
    J.insert( 6, 7 ) = -t34;
    J.insert( 6, 8 ) = -t34;
    J.insert( 6, 9 ) = -t34;
    J.insert( 7, 0 ) = -t45;
    J.insert( 7, 1 ) = -t45;
    J.insert( 7, 2 ) = -1799.0 / 100000000.0 * x( 0 ) - 1799.0 / 100000000.0 * x( 1 ) - 1799.0 / 50000000.0 * x( 2 ) -
                       1799.0 / 100000000.0 * x( 3 ) - 1799.0 / 100000000.0 * x( 4 ) - 1799.0 / 100000000.0 * x( 5 ) -
                       1799.0 / 100000000.0 * x( 6 ) - 1799.0 / 100000000.0 * x( 7 ) - 1799.0 / 100000000.0 * x( 8 ) -
                       1799.0 / 100000000.0 * x( 9 );
    J.insert( 7, 3 ) = x( 7 ) - t45;
    J.insert( 7, 4 ) = -t45;
    J.insert( 7, 5 ) = -t45;
    J.insert( 7, 6 ) = -t45;
    J.insert( 7, 7 ) = x( 3 ) - t45;
    J.insert( 7, 8 ) = -t45;
    J.insert( 7, 9 ) = -t45;
    J.insert( 8, 0 ) = -431.0 / 2000000.0 * t60 - t65;
    J.insert( 8, 1 ) = -t65;
    J.insert( 8, 2 ) = -431.0 / 4000000.0 * t63 *
                       ( x( 0 ) + x( 1 ) + 2.0 * x( 2 ) + x( 3 ) + x( 4 ) + x( 5 ) + x( 6 ) + x( 7 ) + x( 8 ) +
                         x( 9 ) );
    J.insert( 8, 3 ) = x( 8 ) - t65;
    J.insert( 8, 4 ) = -t65;
    J.insert( 8, 5 ) = -t65;
    J.insert( 8, 6 ) = -t65;
    J.insert( 8, 7 ) = -t65;
    J.insert( 8, 8 ) = x( 3 ) - t65;
    J.insert( 8, 9 ) = -t65;
    J.insert( 9, 0 ) = -t74;
    J.insert( 9, 1 ) = -t74;
    J.insert( 9, 2 ) = -t74;
    J.insert( 9, 3 ) = 2.0 * x( 9 ) * x( 3 ) - 1923.0 / 25000000.0 * t15 - t74;
    J.insert( 9, 4 ) = -t74;
    J.insert( 9, 5 ) = -t74;
    J.insert( 9, 6 ) = -t74;
    J.insert( 9, 7 ) = -t74;
    J.insert( 9, 8 ) = -t74;
    J.insert( 9, 9 ) = 49998077.0 / 50000000.0 * t73;

    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 2.998999999984848914379136, 3.988799803896425004703170, 5.999871650821869274792302,
      0.001000000000000000000000000, 0.0002566983540637031238314806, 0.0002969410500484113283243128,
      0.02159005444686014319849639, 1.411134475268817487496289, 5.723908589993461301590475, -7.072216461786208715286129;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 1 );
  }
};
