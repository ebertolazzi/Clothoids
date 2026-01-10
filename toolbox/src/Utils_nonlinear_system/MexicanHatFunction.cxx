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

static inline string ini_msg_MexicanHatFunction( real_type tau )
{
  return fmt::format( "Mexican Hat Function, tau = {}", tau );
}

class MexicanHatFunction : public NonlinearSystem
{
  real_type tau;

public:
  MexicanHatFunction( real_type tau_in )
    : NonlinearSystem(
        ini_msg_MexicanHatFunction( tau_in ),
        "@article{Grippo:1991,\n"
        "  author  = {Grippo, L. and Lampariello, F. and Lucidi, S.},\n"
        "  title   = {A Class of Nonmonotone Stabilization Methods\n"
        "             in Unconstrained Optimization},\n"
        "  journal = {Numer. Math.},\n"
        "  year    = {1991},\n"
        "  volume  = {59},\n"
        "  number  = {1},\n"
        "  pages   = {779--805},\n"
        "  doi     = {10.1007/BF01385810},\n"
        "}\n",
        2 )
    , tau( tau_in )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type t1 = x( 0 ) * x( 0 );
    real_type t2 = x( 1 ) - t1;
    real_type t3 = t2 * t2;
    real_type t6 = power2( 1.0 - x( 0 ) );
    real_type t7 = 10000.0 * t3 + t6 - 0.2E-1;
    f( 0 )       = 2 * ( 1 - x( 0 ) ) + tau * 4.0 * t7 * ( ( 1 - 20000.0 * t2 ) * x( 0 ) - 1 );
    f( 1 )       = 2 * ( 1 - x( 1 ) ) + tau * 40000.0 * t7 * t2;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    // Alias per chiarezza
    real_type const x0 = x( 0 );
    real_type const x1 = x( 1 );

    // Calcoli intermedi con opportuno ordinamento
    real_type const x0_sq        = x0 * x0;
    real_type const one_minus_x0 = 1.0 - x0;
    real_type const t2           = x1 - x0_sq;  // x1 - x0²

    // Calcoli di t7: 10000*(x1-x0²)² + (1-x0)² - 0.02
    real_type const t2_sq           = t2 * t2;
    real_type const one_minus_x0_sq = one_minus_x0 * one_minus_x0;
    real_type const t7              = 10000.0 * t2_sq + one_minus_x0_sq - 0.02;

    // g0 = (1 - 20000*t2)*x0 - 1 = x0 - 20000*x0*t2 - 1
    real_type const g0 = ( 1.0 - 20000.0 * t2 ) * x0 - 1.0;

    // Derivate di t7
    real_type const dt7_dx0 = -40000.0 * t2 * x0 - 2.0 * one_minus_x0;
    real_type const dt7_dx1 = 20000.0 * t2;

    // Derivate di g0 - calcolate in modo stabile
    // dg0/dx0 = 1 - 20000*t2 + 40000*x0²
    real_type const dg0_dx0 = 1.0 - 20000.0 * t2 + 40000.0 * x0_sq;

    // dg0/dx1 = -20000*x0
    real_type const dg0_dx1 = -20000.0 * x0;

    // Termini comuni
    real_type const tau_4     = 4.0 * tau;
    real_type const tau_40000 = 40000.0 * tau;

    // Calcolo dei prodotti in modo stabile (evitando cancellazioni)
    real_type const term1_00 = dt7_dx0 * g0;
    real_type const term2_00 = t7 * dg0_dx0;
    real_type const J00_term = term1_00 + term2_00;

    real_type const term1_01 = dt7_dx1 * g0;
    real_type const term2_01 = t7 * dg0_dx1;
    real_type const J01_term = term1_01 + term2_01;

    // Derivate di t2
    real_type const dt2_dx0 = -2.0 * x0;
    real_type const dt2_dx1 = 1.0;

    real_type const term1_10 = dt7_dx0 * t2;
    real_type const term2_10 = t7 * dt2_dx0;
    real_type const J10_term = term1_10 + term2_10;

    real_type const term1_11 = dt7_dx1 * t2;
    real_type const term2_11 = t7 * dt2_dx1;
    real_type const J11_term = term1_11 + term2_11;

    // Inizializzazione della matrice Jacobiana
    J.resize( n, n );
    J.setZero();

    // Elementi dello Jacobiano
    // Aggiungiamo -2 direttamente agli elementi diagonali
    J.insert( 0, 0 ) = -2.0 + tau_4 * J00_term;
    J.insert( 0, 1 ) = tau_4 * J01_term;
    J.insert( 1, 0 ) = tau_40000 * J10_term;
    J.insert( 1, 1 ) = -2.0 + tau_40000 * J11_term;

    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 3 );
    auto & x0{ x_vec[0] };
    auto & x1{ x_vec[1] };
    auto & x2{ x_vec[2] };
    x0.resize( n );
    x1.resize( n );
    x2.resize( n );
    x0 << 0.86, 0.72;
    x1 << 0.85858, 0.7371534;
    x2 << 1.1414204, 1.3028457;
  }
};
