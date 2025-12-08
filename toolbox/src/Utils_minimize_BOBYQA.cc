/*
 *
 * Implementation of Mike Powell's BOBYQA algorithm for minimizing a function
 * of many variables.  The method is "derivatives free" (only the function
 * values are needed) and accounts for bound constraints on the variables.  The
 * algorithm is described in:
 *
 *   M.J.D. Powell, "The BOBYQA Algorithm for Bound Constrained Optimization
 *   Without Derivatives."  Technical report, Department of Applied Mathematics
 *   and Theoretical Physics, University of Cambridge (2009).
 *
 * The present code is based on the original FORTRAN version written by Mike
 * Powell who kindly provides his code on demand (at mjdp@cam.ac.uk) and has
 * been converted to C by É. Thiébaut.
 *
 * Copyright (c) 2009, Mike Powell (FORTRAN version).
 * Copyright (c) 2015, Éric Thiébaut (C version).
 *
 * Read the accompanying `LICENSE` file for details.
 */

#include <stdio.h>
#include <math.h>

#include "Utils_fmt.hh"
#include "Utils_minimize_BOBYQA.hh"

#define OUTPUT stdout

namespace Utils
{

  using std::abs;
  using std::max;
  using std::min;
  using std::sqrt;

  template <typename Scalar>
  std::string
  BOBYQA_minimizer<Scalar>::print_vec( Vector const & x, integer max_elem ) const
  {
    integer x_size = x.size();

    if ( x_size == 0 ) return "[]";

    std::string result = "[";
    if ( x_size <= max_elem )
    {
      for ( integer i = 0; i < x_size; ++i )
      {
        if ( i != 0 ) result += ", ";
        result += fmt::format( "{:.6}", x[i] );
      }
    }
    else
    {
      for ( int i = 0; i < max_elem; ++i )
      {
        if ( i != 0 ) result += ", ";
        result += fmt::format( "{}", x[i] );
      }
      result += ", ...";
      // Nota: non mostriamo gli elementi finali in questo caso, ma potremmo se volessimo
      // l'esempio fornito invece mostra solo i primi max_elem e poi "..."
    }
    result += "]";
    return result;
  }

  template <typename Scalar>
  Scalar
  BOBYQA_minimizer<Scalar>::eval( Vector const & x ) {
    ++m_num_f_eval;
    // Controlla se sono state esaurite le valutazioni durante il rescue
    if ( m_num_f_eval > m_max_f_eval ) {
      m_status = Status::BOBYQA_TOO_MANY_EVALUATIONS;
      m_reason = fmt::format( "BOBYQA_minimizer::eval( x ) has been called MAXFUN={} times", m_max_f_eval );
      UTILS_ERROR( m_reason );
    }
    Scalar f = m_fun( x );
    if ( m_print_level == 3 )
      fmt::print( "    Function n.{} F(X) = {:15}   X: {}\n", m_num_f_eval, fmt::format( "{:.9}", f ),
                  print_vec( x, 6 ) );
    return f;
  }

  /**
   * @brief Hessian-vector product for the model (linear-algebra view).
   *
   * For the quadratic model used by BOBYQA the Hessian has two parts:
   *   H = HQ (explicit symmetric NxN matrix) + sum_{k=1..NPT} pq[k] * xpt_k xpt_k^T
   *
   * This routine computes hs = H * s for a given vector s by:
   *   1) hs := HQ * s        (HQ is stored as a symmetric dense matrix m_HQ)
   *   2) hs += m_xpt * ( m_pq.asDiagonal() * ( m_xpt^T * s ) )
   *
   * Linear algebra identities:
   *   let y = m_xpt^T * s   (length NPT)
   *   then second term = m_xpt * ( diag(pq) * y ) = sum_k pq[k] * y[k] * xpt_k
   *
   * Fortran mapping: corresponds to contributions computed from HQ and PQ/XPT
   * in subroutines TRSBOX, UPDATE and other places where H*s is needed.
   */

  template <typename Scalar>
  void
  BOBYQA_minimizer<Scalar>::compute_hessian_product( Vector const & s, Vector & hs ) const
  {
    // 2. Parte XPT/PQ
    // pb: PQ has NPT values, m_xpt is N × NPT (trasposta del Fortran)
    auto HQ = m_HQ.template selfadjointView<Eigen::Upper>();
    if ( m_npt <= m_num_f_eval ) {
      hs.noalias() = HQ * s + m_xpt * (m_pq.array() * (m_xpt.transpose() * s).array()).matrix();
    } else {
      hs.noalias() = HQ * s;
    }
  }

  template <typename Scalar>
  void
  BOBYQA_minimizer<Scalar>::add_hessian_product( Vector const & s, Vector & hs ) const
  {
    // 2. Parte XPT/PQ
    // pb: PQ has NPT values, m_xpt is N × NPT (trasposta del Fortran)
    auto HQ = m_HQ.template selfadjointView<Eigen::Upper>();
    if ( m_npt <= m_num_f_eval ) {
      hs.noalias() += HQ * s + m_xpt * (m_pq.array() * (m_xpt.transpose() * s).array()).matrix();
    } else {
      hs.noalias() += HQ * s;
    }
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Main minimization routine for BOBYQA algorithm.
   *
   * This function initializes the BOBYQA minimizer with problem dimensions,
   * validates parameters, processes bound constraints, and calls the core
   * optimization routine.
   *
   * @tparam Scalar Floating point type (e.g., double, float)
   * @param N Dimension of the problem (number of variables)
   * @param NPT Number of interpolation points (must satisfy N+2 ≤ NPT ≤ (N+1)(N+2)/2)
   * @param objfun Objective function to minimize
   * @param[in,out] X On input: initial guess. On output: final solution.
   * @param XL Vector of lower bounds
   * @param XU Vector of upper bounds
   * @return Status::BOBYQA_SUCCESS on success, appropriate error code otherwise
   *
   * @throws No exceptions thrown, but returns error status codes
   *
   * @details
   * The function performs the following steps:
   * 1. Stores problem dimensions and bounds
   * 2. Resizes all internal vectors and matrices
   * 3. Validates NPT parameter
   * 4. Checks bound feasibility
   * 5. Adjusts initial point if too close to bounds
   * 6. Calls bobyqb() for core optimization
   */
  template <typename Scalar>
  typename BOBYQA_minimizer<Scalar>::Status
  BOBYQA_minimizer<Scalar>::minimize( integer const         N,
                                      integer const         NPT,
                                      bobyqa_objfun const & fun,
                                      Vector &              X,
                                      Vector const &        XL,
                                      Vector const &        XU )
  {
    // Call the core BOBYQB optimization routine with adjusted starting point

    using Eigen::Index;
    using std::max;
    using std::min;

    m_fun = fun;

    // Store problem dimensions
    m_nv   = N;
    m_npt  = NPT;
    m_dim  = m_npt + m_nv;
    m_nptm = m_npt - m_nv - 1;  // Numero colonne di ZMAT (NPTM)

    // Validate NPT parameter (as in original BOBYQA implementation)
    const integer np = m_nv + 1;
    if ( m_npt < m_nv + 2 || m_npt > ( m_nv + 2 ) * np / 2 )
    {
      print_error( "NPT is not in the required interval" );
      return Status::BOBYQA_BAD_NPT;
    }

    // Store bounds (Eigen uses copy-on-write for assignments)
    m_x_lower = XL;
    m_x_upper = XU;

    // Resize all internal storage using Eigen
    // Note: All resizing operations are O(1) for fixed-size vectors,
    // and only allocate memory when size changes for dynamic vectors
    m_xbdi.resize( m_nv );
    m_x_base.resize( m_nv );
    m_x_new.resize( m_nv );
    m_g_new.resize( m_nv );
    m_x_alt.resize( m_nv );
    m_f_val.resize( m_npt );
    m_x_opt.resize( m_nv );
    m_g_opt.resize( m_nv );
    m_pq.resize( m_npt );
    m_s_lower.resize( m_nv );
    m_s_upper.resize( m_nv );
    m_d.resize( m_nv );
    m_v_lag.resize( m_dim );
    m_g_lag.resize( m_nv );  // Used with size m_nv in other routines

    m_hcol.resize( m_npt );
    m_curv.resize( m_nv );

    m_HQ.resize( m_nv, m_nv );
    m_xpt.resize( m_nv, m_npt );  // Matrix: N rows, NPT columns
    m_B.resize( m_nv, m_dim );    // Matrix: N rows, (NPT+N) columns

    // ZMAT matrix: NPT rows, (NPT-N-1) columns (minimum 1 column)
    m_Z.resize( m_npt, m_nptm );

    m_ptsaux.resize( 2, m_nv );  // 2 rows, N columns
    m_ptsid.resize( m_npt );

    m_s.resize( m_nv );
    m_hs.resize( m_nv );
    m_hred.resize( m_nv );

    // working array
    m_WNPT.resize( m_npt );

    // Compute bound differences and check feasibility using vectorized operations
    Scalar const min_diff = (m_x_upper - m_x_lower).minCoeff();  // XU - XL

    // Check if any bound difference is less than 2*RHOBEG
    if ( min_diff < Scalar( 2 ) * m_rhobeg )
    {
      print_error( "one of the differences XU(I)-XL(I) is less than 2*RHOBEG" );
      return Status::BOBYQA_TOO_CLOSE;
    }

    // Create per-component slack variables and adjust X according to the Fortran logic
    for ( integer j = 0; j < m_nv; ++j )
    {
      auto const & L  = m_x_lower(j);
      auto const & U  = m_x_upper(j);
      auto       & Xj = X( j );

      Scalar temp = U - L;
      Scalar wsl  = L - Xj;  // xl - x
      Scalar wsu  = U - Xj;  // xu - x

      if ( wsl >= -m_rhobeg )
      {
        if ( wsl >= Scalar( 0 ) )
        {
          Xj  = L;
          wsl = Scalar( 0 );
          wsu = temp;
        }
        else
        {
          Xj  = L + m_rhobeg;
          wsl = -m_rhobeg;
          wsu = std::max( U - Xj, m_rhobeg );
        }
      }
      else if ( wsu <= m_rhobeg )
      {
        if ( wsu <= Scalar( 0 ) )
        {
          X( j ) = U;
          wsl    = -temp;
          wsu    = Scalar( 0 );
        }
        else
        {
          X( j ) = U - m_rhobeg;
          wsl    = std::min( L - Xj, -m_rhobeg );
          wsu    = m_rhobeg;
        }
      }

      m_s_lower( j ) = wsl;
      m_s_upper( j ) = wsu;
    }

    // The initial X may have changed; ensure it lies within bounds
    X = X.cwiseMax( m_x_lower ).cwiseMin( m_x_upper );

    // workspace
    //Vector W1 = Vector::Zero( std::max( 5 * m_nv, 2 * m_npt ) );

    // many scalars used by algorithm
    Scalar ratio = 0;

    m_status = Status::BOBYQA_SUCCESS;

    // INITIALIZATION (calls prelim which must fill m_xpt, m_f_val, etc.)
    prelim( X );

    // If prelim didn't have enough function evaluations then abort like Fortran
    if ( m_num_f_eval < m_npt )
    {
      if ( m_print_level > 0 ) print_error( "CALFUN has been called MAXFUN times" );
      return Status::BOBYQA_TOO_MANY_EVALUATIONS;
    }

    // set m_x_opt from m_xpt row m_kopt (0-based)
    m_x_opt               = m_xpt.col( m_kopt );
    Scalar m_x_opt_square = m_x_opt.squaredNorm();

    Scalar  fsave  = m_f_val( 0 );
    integer kbase  = 0;
    integer ntrits = 0;
    Scalar  diffa  = 0;
    Scalar  diffb  = 0;
    Scalar  diffc  = 0;
    Scalar  f      = 0;
    Scalar  dnorm  = 0;
    integer itest  = 0;

    m_rho    = m_rhobeg;
    m_crvmin = 0;
    m_dsq    = 0;
    m_alpha  = 0;
    m_beta   = 0;
    m_delta  = m_rhobeg;
    m_cauchy = 0;
    m_denom  = 0;
    m_adelt  = 0;
    m_knew   = -1;
    m_distsq = 0;

    n_num_f_saved  = m_num_f_eval;
    n_num_f_rescue = m_num_f_eval;
    
    m_s.setZero();

    // Main state-machine loop variables: we'll use an enum and switch
    enum class Phase
    {
      UPDATE_GRADIENT,
      TRUST_REGION,
      SHIFT_BASE,
      RESCUE,
      ALTMOV,
      COMPUTE_VLAG,
      EVALUATE,
      FIND_FAR,
      REDUCE_RHO,
      DONE,
      ERROR
    };

    auto info = [&]( string const & name ) -> void {
      return;
      std::cout << "\n\n" << name << "\n";
      std::cout << "m_num_f_eval   = " << m_num_f_eval << '\n';
      std::cout << "m_knew         = " << m_knew << '\n';
      std::cout << "m_rho          = " << m_rho << '\n';
      std::cout << "m_crvmin       = " << m_crvmin << '\n';
      std::cout << "m_dsq          = " << m_dsq << '\n';
      std::cout << "m_alpha        = " << m_alpha << '\n';
      std::cout << "m_beta         = " << m_beta << '\n';
      std::cout << "m_delta        = " << m_delta << '\n';
      std::cout << "m_adelt        = " << m_adelt << '\n';
      std::cout << "m_cauchy       = " << m_cauchy << '\n';
      std::cout << "m_denom        = " << m_denom << '\n';
      std::cout << "m_x_opt_square = " << m_x_opt_square << '\n';
      std::cout << "m_distsq = " << m_distsq << '\n';
      std::cout << "m_B =\n"     << m_B << '\n';
      std::cout << "m_Z =\n"     << m_Z << '\n';
      std::cout << "m_x_new =\n" << m_x_new.transpose() << '\n';
      std::cout << "m_g_new =\n" << m_g_new.transpose() << '\n';
    };

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Lambda per l'aggiornamento del gradiente nel punto ottimo corrente.
     *
     * Calcola il gradiente quadratico m_g_opt nel punto ottimo m_x_opt,
     * includendo sia la parte quadratica (m_hq) che la parte polinomiale (m_pq).
     * Viene eseguita quando il punto ottimo cambia rispetto alla base corrente.
     */
    auto phase_update_gradient = [&]() -> Phase
    {
      if ( m_kopt != kbase ) add_hessian_product( m_x_opt, m_g_opt );
      return Phase::TRUST_REGION;  // Passa alla fase successiva
    };

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Lambda per la risoluzione del problema della regione di fiducia.
     *
     * Risolve il problema quadratico vincolato nella regione di fiducia usando trsbox().
     * Determina se il passo è sufficientemente grande o se è necessario ridurre rho
     * o cercare un punto di interpolazione lontano.
     */
    auto phase_trust_region = [&]() -> Phase
    {
      // Risolve il problema quadratico vincolato nella regione di fiducia
      trsbox();

      // Calcola norma del passo e la limita alla dimensione della regione di fiducia
      dnorm = min( sqrt( m_dsq ), m_delta );

      // Controlla se il passo è troppo piccolo rispetto a rho
      if ( dnorm < Scalar( 0.5 ) * m_rho )
      {
        ntrits   = -1;                               // Marca che non è stata fatta una valutazione di funzione
        m_distsq = Scalar( 100 ) * power2( m_rho );  // Soglia per cercare punto lontano

        // Controlli per decidere se ridurre rho o cercare punto lontano
        if ( m_num_f_eval <= n_num_f_saved + 2 ) return Phase::FIND_FAR;

        Scalar errbig = std::max( { diffa, diffb, diffc } );  // Errore di approssimazione
        Scalar frhosq = Scalar( 0.125 ) * power2( m_rho );

        // Se l'errore è grande rispetto alla curvatura minima, cerca punto lontano
        if ( m_crvmin > 0 && errbig > frhosq * m_crvmin ) return Phase::FIND_FAR;

        // Verifica condizioni per ridurre ρ basate sulle condizioni di ottimalità
        Scalar const bdtol    = errbig / m_rho;
        Scalar const half_rho = Scalar( 0.5 ) * m_rho;
        bool         reduce   = true;

        // Precomputo curvatura diagonale HQ + Σ pq * xpt^2 (Eigen-optimizzato)
        m_curv = m_HQ.diagonal();  // m_nv
        m_curv.noalias() += ( m_xpt.array().square().rowwise() * m_pq.transpose().array() ).rowwise().sum().matrix();
        for ( integer j = 0; j < m_nv; ++j )
        {
          auto const & xnj = m_x_new( j );
          auto const & L   = m_s_lower( j );
          auto const & U   = m_s_upper( j );
          Scalar bdtest = bdtol;
          // Direzione bloccata sul bound → use slack distances like Fortran
          if ( xnj == L ) bdtest = L;
          if ( xnj == U ) bdtest = -U;

          // If slack indicates potential improvement then adjust by curvature
          if ( bdtest < bdtol )
          {
            // Aggiungi metà curvatura * rho
            bdtest += half_rho * m_curv( j );

            if ( bdtest < bdtol )
            {
              reduce = false;  // Non ridurre rho → serve punto lontano
              break;
            }
          }
        }
        return reduce ? Phase::REDUCE_RHO : Phase::FIND_FAR;
      }

      ++ntrits;                  // Incrementa contatore valutazioni nella regione di fiducia
      return Phase::SHIFT_BASE;  // Passa alla fase di shift della base
    };

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Lambda per lo shift della base e l'aggiornamento delle matrici.
     *
     * Quando il passo è piccolo rispetto alla distanza dal centro, sposta la base
     * al punto ottimo corrente e aggiorna tutte le matrici (BMAT, ZMAT, HQ)
     * per mantenere l'interpolazione quadratica.
     */
    auto phase_shift_base = [&]() -> Phase
    {
      // Controlla se il passo è significativamente più piccolo della distanza dal centro
      if ( m_dsq <= m_x_opt_square * Scalar( 0.001 ) )
      {
        Scalar fracsq = m_x_opt_square * Scalar( 0.25 );
        Scalar sumpq  = m_pq.sum();

        Vector wn( m_nv );

        // Fase 1: Aggiornamento di BMAT con contributi dai punti di interpolazione
        for ( integer k = 0; k < m_npt; ++k )
        {
          Scalar sum  = m_xpt.col( k ).dot( m_x_opt ) - Scalar( 0.5 ) * m_x_opt_square;
          m_WNPT( k ) = sum;
          Scalar temp = fracsq - Scalar( 0.5 ) * sum;

          // salva colonna in temporaneo
          wn.noalias() = m_B.col( k );

          // m_v_lag usato come temporaneo
          m_v_lag.head( m_nv ) = sum * m_xpt.col( k ) + temp * m_x_opt;

          // Aggiorna la parte inferiore di BMAT
          for ( integer i = 0; i < m_nv; ++i )
            m_B.col( m_npt + i ).head( i + 1 ) += wn( i ) * m_v_lag.head( i + 1 ) + m_v_lag( i ) * wn.head( i + 1 );
        }

        // Fase 2
        for ( integer jj = 0; jj < m_nptm; ++jj )
        {
          m_v_lag.head( m_npt ) = m_WNPT.cwiseProduct( m_Z.col( jj ) );
          Scalar sumw           = m_v_lag.head( m_npt ).sum();
          Scalar sumz           = m_Z.col( jj ).sum();

          for ( integer j = 0; j < m_nv; ++j )
          {
            Scalar sum = ( fracsq * sumz - Scalar( 0.5 ) * sumw ) * m_x_opt( j );
            sum += m_v_lag.head( m_npt ).dot( m_xpt.row( j ) );
            wn( j ) = sum;
            m_B.row( j ).head( m_npt ) += sum * m_Z.col( jj );
          }

          for ( integer i = 0; i < m_nv; ++i ) m_B.col( i + m_npt ).head( i + 1 ) += wn( i ) * wn.head( i + 1 );
        }

        // wn.head(n) = m_xpt * m_pq - m_half * sumpq * m_x_opt
        wn = m_xpt * m_pq - ( Scalar( 0.5 ) * sumpq ) * m_x_opt;

        // m_xpt = m_xpt - spread(m_x_opt, 2, m_npt)
        m_xpt.colwise() -= m_x_opt;

        // Fase 3
        for ( integer j = 0; j < m_nv; ++j )
        {
          auto & xoj = m_x_opt( j );
          auto & Wj  = wn( j );

          for ( integer i = 0; i <= j; ++i )
          {
            m_HQ( i, j ) += wn( i ) * xoj + m_x_opt( i ) * Wj;
            // Rende simmetrica la parte di BMAT
            m_B( j, m_npt + i ) = m_B( i, m_npt + j );
          }
        }

        // Fase 4: Shift completo della base e aggiornamento dei bound
        m_x_base += m_x_opt;   // Sposta la base
        m_x_new -= m_x_opt;    // Corregge il nuovo punto
        m_s_lower -= m_x_opt;  // Aggiorna limite inferiore relativo
        m_s_upper -= m_x_opt;  // Aggiorna limite superiore relativo
        m_x_opt.setZero();     // Reset del punto ottimo relativo
        m_x_opt_square = 0;    // Reset della norma al quadrato
      }

      return ( ntrits == 0 ) ? Phase::ALTMOV : Phase::COMPUTE_VLAG;
    };

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Lambda per l'operazione di rescue.
     *
     * Chiama la procedura rescue() quando si sospettano problemi numerici
     * o quando i denominatori diventano troppo piccoli. Rigenera il set
     * di punti di interpolazione mantenendo l'interpolazione quadratica.
     */
    auto phase_rescue = [&]() -> Phase
    {
      n_num_f_saved = m_num_f_eval;  // Salva il numero corrente di valutazioni
      kbase         = m_kopt;        // Salva il punto ottimo corrente come base

      // Esegue la procedura di rescue
      rescue();
      m_x_opt_square = 0;

      // Se il punto ottimo è cambiato, aggiorna m_x_opt e m_x_opt_square
      if ( m_kopt != kbase )
      {
        m_x_opt        = m_xpt.col( m_kopt );
        m_x_opt_square = m_x_opt.squaredNorm();
      }

      // Controlla se sono state esaurite le valutazioni durante il rescue
      if ( m_num_f_eval >= m_max_f_eval )
      {
        m_num_f_eval = m_max_f_eval;
        m_reason     = "CALFUN has been called MAXFUN times";
        return Phase::ERROR;
      }

      n_num_f_rescue = m_num_f_eval;  // Aggiorna contatore rescue

      // Decide la fase successiva in base allo stato
      if ( n_num_f_saved < m_num_f_eval )
      {
        n_num_f_saved = m_num_f_eval;
        return Phase::UPDATE_GRADIENT;  // Nuove valutazioni, aggiorna gradiente
      }
      else if ( ntrits > 0 )
      {
        return Phase::TRUST_REGION;  // Ritorna alla regione di fiducia
      }
      else
      {
        return Phase::ALTMOV;  // Prova movimento alternativo
      }
    };

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Lambda per il calcolo del movimento alternativo.
     *
     * Chiama altmov() per generare un passo alternativo quando la regione
     * di fiducia è ristretta o quando non sono state fatte valutazioni di funzione.
     */
    auto phase_altmov = [&]() -> Phase
    {
      // Calcola un passo alternativo
      altmov();

      // Calcola la direzione d dal punto ottimo al nuovo punto
      m_d = m_x_new - m_x_opt;

      return Phase::COMPUTE_VLAG;  // Passa al calcolo dei coefficienti di Lagrange
    };

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Lambda per il calcolo dei coefficienti di Lagrange.
     *
     * Calcola i coefficienti di Lagrange m_v_lag per il punto di interpolazione
     * da sostituire e il denominatore m_denom per la regola di aggiornamento.
     * Sceglie anche il punto da rimuovere (m_knew) se necessario.
     */
    auto phase_compute_vlag = [&]() -> Phase
    {
      Vector W0( m_npt );
      // Fase 1: Calcola contributi base dai punti di interpolazione
      for ( integer k = 0; k < m_npt; ++k )
      {
        Scalar suma    = m_xpt.col( k ).dot( m_d );               // Proiezione sulla direzione d
        Scalar sumb    = m_xpt.col( k ).dot( m_x_opt );           // Proiezione sul punto ottimo
        Scalar sum     = m_B.col( k ).dot( m_d );                 // Contributo di BMAT
        W0( k )        = suma * ( Scalar( 0.5 ) * suma + sumb );  // Termine quadratico
        m_v_lag( k )   = sum;                                     // Coefficiente di Lagrange
        m_WNPT( k )    = suma;                                    // Salva per uso successivo
      }

      // Fase 2: Calcola contributo da ZMAT (parte ortogonale)
      m_beta = 0;
      for ( integer jj = 0; jj < m_nptm; ++jj )
      {
        Scalar sum = m_Z.col( jj ).dot( W0 );
        m_beta -= sum * sum;                           // Aggiorna beta (contributo negativo)
        m_v_lag.head( m_npt ) += sum * m_Z.col( jj );  // Aggiorna coefficienti di Lagrange
      }

      // Fase 3: Calcola norme e prodotti scalari
      m_dsq = m_d.squaredNorm();  // Norma al quadrato della direzione

      // dx = dot(m_d, m_x_opt)
      Scalar dx = m_d.dot( m_x_opt );

      // m_v_lag(m_npt : m_npt+n-1) = m_B(0:n-1, 0:m_npt-1) * wn
      m_v_lag.segment( m_npt, m_nv ) = m_B.topLeftCorner( m_nv, m_npt ) * W0;

      // bsum = dot(m_d, m_v_lag_tail)
      Scalar bsum = m_d.dot( m_v_lag.segment( m_npt, m_nv ) );

      // m_v_lag_tail += transpose(m_B(:, m_npt:m_npt+n-1)) * m_d
      m_v_lag.segment( m_npt, m_nv ) += m_B.topRightCorner( m_nv, m_nv ).transpose() * m_d;

      // bsum += dot(m_d, m_v_lag_tail)
      bsum += m_d.dot( m_v_lag.segment( m_npt, m_nv ) );

      // Calcola beta finale (distanza quadratica normalizzata)
      m_beta = dx * dx + m_dsq * ( m_x_opt_square + 2 * dx + Scalar( 0.5 ) * m_dsq ) + m_beta - bsum;
      m_v_lag( m_kopt ) += 1;  // Condizione di interpolazione nel punto ottimo

      // Fase 4: Scelta del denominatore e controllo di validità
      if ( ntrits == 0 )
      {
        // Prima iterazione nella regione di fiducia
        m_denom = power2( m_v_lag( m_knew ) ) + m_alpha * m_beta;

        // Controlla se il denominatore è peggiore del passo di Cauchy
        if ( m_denom < m_cauchy && m_cauchy > 0 )
        {
          // Revert al passo alternativo
          m_x_new  = m_x_alt;
          m_d      = m_x_new - m_x_opt;
          m_cauchy = 0;
          return Phase::COMPUTE_VLAG;  // Ricalcola con passo alternativo
        }

        // Controlla cancellazione numerica nel denominatore
        if ( m_denom <= Scalar( 0.5 ) * power2( m_v_lag( m_knew ) ) )
        {
          if ( m_num_f_eval > n_num_f_rescue ) return Phase::RESCUE;  // Richiede rescue per problemi numerici
          m_reason = "of much cancellation in a denominator";
          return Phase::ERROR;  // Errore irreparabile
        }
      }
      else
      {
        // Iterazioni successive: sceglie il punto da rimuovere (m_knew)
        Scalar delsq  = power2( m_delta );
        Scalar scaden = 0;  // Denominatore scalato massimo
        Scalar biglsq = 0;  // Termine di confronto
        // m_knew        = 0; non devo azzerare se poi m_knew non viene aggiornato

        for ( integer k = 0; k < m_npt; ++k )
        {
          if ( k == m_kopt ) continue;  // Salta il punto ottimo

          // Calcola hdiag (contributo ortogonale)
          Scalar hdiag = m_Z.row( k ).squaredNorm();
          // for ( integer jj = 0; jj < m_nptm; ++jj ) { hdiag += m_zmat( k, jj ) * m_zmat( k, jj ); }

          // Calcola denominatore per questo punto
          Scalar den = m_beta * hdiag + m_v_lag( k ) * m_v_lag( k );

          // Calcola distanza normalizzata dal punto ottimo
          m_distsq    = ( m_xpt.col( k ) - m_x_opt ).squaredNorm();
          Scalar temp = m_distsq / delsq;
          temp       *= temp;
          if ( temp < 1 ) temp = 1;  // Evita valori troppo piccoli

          // Aggiorna scaden e biglsq
          if ( temp * den > scaden )
          {
            scaden  = temp * den;
            m_knew  = k;  // 0-based index
            m_denom = den;
          }
          temp *= power2( m_v_lag( k ) );
          biglsq = max( biglsq, temp );
        }

        // Controlla se il denominatore è sufficientemente grande
        if ( scaden <= Scalar( 0.5 ) * biglsq )
        {
          if ( m_num_f_eval > n_num_f_rescue ) return Phase::RESCUE;  // Problemi numerici, richiede rescue
          m_reason = "of much cancellation in a denominator";
          return Phase::ERROR;  // Errore irreparabile
        }
      }

      return Phase::EVALUATE;  // Procedi con la valutazione della funzione
    };

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Lambda per la valutazione della funzione obiettivo.
     *
     * Valuta la funzione obiettivo nel nuovo punto, aggiorna il modello quadratico,
     * e decide se accettare il passo, ridurre la regione di fiducia, o ridurre rho.
     */
    auto phase_evaluate = [&]() -> Phase
    {
      // Fase 1: Proietta il nuovo punto sui bound con trattamento speciale
      X = ( m_x_base + m_x_new ).cwiseMax( m_x_lower ).cwiseMin( m_x_upper );

      // Controllo limite massimo valutazioni
      if ( m_num_f_eval >= m_max_f_eval )
      {
        m_reason = "CALFUN has been called MAXFUN times";
        return Phase::ERROR;
      }

      // Valutazione funzione obiettivo
      f = eval( X );

      // Caso speciale: nessuna valutazione di funzione nella regione di fiducia
      if ( ntrits == -1 )
      {
        fsave = f;
        return Phase::DONE;
      }

      // Fase 2: Calcolo riduzione predetta dal modello quadratico
      Scalar fopt  = m_f_val( m_kopt );
      Scalar vquad = m_d.dot( m_g_opt );  // Riduzione predetta: Contributo lineare

      for ( integer j = 0; j < m_nv; ++j )
      {
        for ( integer i = 0; i <= j; ++i )
        {
          Scalar temp = m_d( i ) * m_d( j );
          if ( i == j ) temp = Scalar( 0.5 ) * temp;
          vquad += m_HQ( i, j ) * temp;  // Contributo quadratico
        }
      }

      vquad += Scalar( 0.5 ) * ( m_pq.array() * m_WNPT.array().square() ).sum();
      // for ( integer k = 0; k < m_npt; ++k )
      //   vquad += Scalar( 0.5 ) * m_pq( k ) * power2( W( m_npt + k ) );

      // Fase 3: Aggiorna statistiche errore approssimazione
      Scalar diff = f - fopt - vquad;
      diffc       = diffb;
      diffb       = diffa;
      diffa       = std::abs( diff );

      if ( dnorm > m_rho ) n_num_f_saved = m_num_f_eval;

      // Fase 4: Gestione risultato valutazione
      if ( ntrits > 0 )
      {
        // Controlla se il modello predice correttamente
        if ( vquad >= 0 )
        {
          m_reason = "a trust region step has failed to reduce Q";
          return Phase::ERROR;
        }

        // Calcola ratio riduzione effettiva/predetta
        ratio = ( f - fopt ) / vquad;

        // Regola dimensione regione di fiducia in base al ratio
        if ( ratio <= Scalar( 0.1 ) )
          m_delta = std::min( Scalar( 0.5 ) * m_delta, dnorm );
        else if ( ratio <= Scalar( 0.7 ) )
          m_delta = std::max( Scalar( 0.5 ) * m_delta, dnorm );
        else
          m_delta = std::max( Scalar( 0.5 ) * m_delta, 2 * dnorm );

        if ( m_delta <= Scalar( 1.5 ) * m_rho ) m_delta = m_rho;

        // Se il passo migliora la funzione, rivaluta scelta punto da rimuovere
        if ( f < fopt )
        {
          integer ksav   = m_knew;
          Scalar  densav = m_denom;
          Scalar  delsq  = m_delta * m_delta;
          Scalar  scaden = 0;
          Scalar  biglsq = 0;
          m_knew         = -1;

          for ( integer k = 0; k < m_npt; ++k )
          {
            Scalar hdiag = m_Z.row( k ).squaredNorm();
            Scalar den = m_beta * hdiag + power2( m_v_lag( k ) );

            m_distsq    = ( m_xpt.col( k ) - m_x_new ).squaredNorm();
            Scalar temp = m_distsq / delsq;
            temp        = temp * temp;
            if ( temp < 1 ) temp = 1;

            if ( scaden < temp * den )
            {
              scaden  = temp * den;
              m_knew  = k;
              m_denom = den;
            }
            biglsq = max( biglsq, temp * power2( m_v_lag( k ) ) );
          }

          if ( scaden <= Scalar( 0.5 ) * biglsq || m_knew == -1 )
          {
            m_knew  = ksav;
            m_denom = densav;
          }
        }
      }

      // Fase 5: Aggiornamento del modello quadratico
      update();  // Aggiorna BMAT, ZMAT

      // Aggiorna HQ con contributo del punto rimosso
      Scalar pqold   = m_pq( m_knew );
      m_pq( m_knew ) = 0;
      for ( integer i = 0; i < m_nv; ++i )
      {
        Scalar temp = pqold * m_xpt( i, m_knew );
        for ( integer j = 0; j <= i; ++j ) { m_HQ( j, i ) += temp * m_xpt( j, m_knew ); }
      }

      // Aggiorna PQ con differenza funzione
      for ( integer jj = 0; jj < m_nptm; ++jj )
      {
        Scalar temp = diff * m_Z( m_knew, jj );
        m_pq += temp * m_Z.col( jj );
      }

      // Incorpora nuovo punto di interpolazione
      m_f_val( m_knew )   = f;
      m_xpt.col( m_knew ) = m_x_new;
      Vector Bkn = m_B.col( m_knew );

      // Aggiorna gradiente con nuovo punto
      for ( integer k = 0; k < m_npt; ++k )
      {
        Scalar suma = m_Z.row( m_knew ).dot( m_Z.row( k ) );
        Scalar sumb = m_xpt.col( k ).dot( m_x_opt );
        Scalar temp = suma * sumb;
        Bkn += temp * m_xpt.col( k );
      }

      m_g_opt += diff * Bkn;

      // Fase 6: Aggiorna punto ottimo se migliorato
      if ( f < fopt )
      {
        m_kopt         = m_knew;
        m_x_opt        = m_x_new;
        m_x_opt_square = m_x_opt.squaredNorm();
        add_hessian_product( m_d, m_g_opt );
      }

      // Fase 7: Controllo qualità interpolazione e aggiornamento forzato
      if ( ntrits > 0 )
      {
        // Calcola interpolante di Frobenius minimo
        Vector wt0( m_npt );
        wt0.setZero();

        for ( integer k = 0; k < m_npt; ++k ) { m_v_lag( k ) = m_f_val( k ) - m_f_val( m_kopt ); }

        for ( integer jj = 0; jj < m_nptm; ++jj )
        {
          Scalar sum = m_Z.col( jj ).dot( m_v_lag.head( m_npt ) );
          wt0 += sum * m_Z.col( jj );
        }

        m_WNPT.noalias() = wt0;
        for ( integer k = 0; k < m_npt; ++k )
        {
          Scalar sum = m_xpt.col( k ).dot( m_x_opt );
          wt0( k ) *= sum;
        }

        // Verifica qualità gradiente interpolante
        Scalar gqsq = 0;  // Norma gradiente quadratico
        Scalar gisq = 0;  // Norma gradiente interpolante

        for ( integer i = 0; i < m_nv; ++i )
        {
          Scalar sum = m_B.row( i ).head( m_npt ).dot( m_v_lag.head( m_npt ) ) + m_xpt.row( i ).dot( wt0 );

          // Gestione vincoli attivi
          auto const & gi = m_g_opt( i );
          auto const & xi = m_x_opt( i );
          if ( xi == m_s_lower( i ) )
          {
            if ( gi < 0 ) gqsq += power2( gi );
            if ( sum < 0 ) gisq += power2( sum );
          }
          else if ( xi == m_s_upper( i ) )
          {
            if ( gi > 0 ) gqsq += power2( gi );
            if ( sum > 0 ) gisq += power2( sum );
          }
          else
          {
            gqsq += power2( gi );
            gisq += power2( sum );
          }

          m_v_lag( m_npt + i ) = sum;
        }

        // Controllo convergenza interpolante
        ++itest;
        if ( gqsq < Scalar( 10 ) * gisq ) itest = 0;

        // Se necessario, sostituisce modello con interpolante minimo
        if ( itest >= 3 )
        {
          m_g_opt = m_v_lag.tail( m_nv );
          m_pq    = m_WNPT;
          m_HQ.setZero();
          itest = 0;
        }
      }

      // Fase 8: Decisione fase successiva
      if ( ntrits == 0 || f <= fopt + Scalar( 0.1 ) * vquad ) return Phase::TRUST_REGION;  // Continua con regione di fiducia
      // Cerca punto lontano per migliorare geometria
      m_distsq = max( power2( 2 * m_delta ), power2( Scalar( 10 ) * m_rho ) );
      return Phase::FIND_FAR;
    };

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Lambda per la ricerca di un punto di interpolazione lontano.
     *
     * Cerca il punto di interpolazione più lontano dal punto ottimo corrente
     * per sostituirlo con un nuovo punto che migliori la geometria dell'interpolazione.
     */
    auto phase_find_far = [&]() -> Phase
    {
      m_knew = -1;

      // Trova il punto più lontano dal punto ottimo
      for ( integer k = 0; k < m_npt; ++k )
      {
        Scalar sum = ( m_xpt.col( k ) - m_x_opt ).squaredNorm();
        if ( sum > m_distsq )
        {
          m_distsq = sum;
          m_knew   = k;  // 0-based index
        }
      }

      if ( m_knew >= 0 )
      {
        // Calcola distanza e adatta parametri
        Scalar dist = sqrt( m_distsq );

        if ( ntrits == -1 )
        {
          // Prima iterazione: riduci delta
          m_delta = min( Scalar( 0.1 ) * m_delta, Scalar( 0.5 ) * dist );
          if ( m_delta <= Scalar( 1.5 ) * m_rho ) m_delta = m_rho;
        }

        ntrits  = 0;  // Reset contatore valutazioni
        m_adelt = max( min( Scalar( 0.1 ) * dist, m_delta ), m_rho );
        m_dsq   = m_adelt * m_adelt;
        return Phase::SHIFT_BASE;  // Procedi con shift base per nuovo punto
      }
      else
      {
        // Nessun punto lontano trovato
        if ( ntrits == -1 ) return Phase::REDUCE_RHO;  // Riduci rho se nessuna valutazione
        if ( ratio > 0 || max( m_delta, dnorm ) > m_rho ) return Phase::TRUST_REGION;  // Continua con regione di fiducia
        return Phase::REDUCE_RHO;  // Riduci rho
      }
    };

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Lambda per la riduzione del parametro rho.
     *
     * Riduce il parametro rho (dimensione caratteristica) quando l'algoritmo
     * ha convergentro alla precisione desiderata o quando è necessario
     * esplorare su scala più fine.
     */
    auto phase_reduce_rho = [&]() -> Phase
    {
      // Controlla se siamo al limite di precisione
      if ( m_rho > m_rhoend )
      {
        // Riduci rho in base al rapporto con rhoend
        m_delta = Scalar( 0.5 ) * m_rho;
        ratio   = m_rho / m_rhoend;

        if ( ratio <= Scalar( 16 ) )
          m_rho = m_rhoend;
        else if ( ratio <= Scalar( 250 ) )
          m_rho = sqrt( ratio ) * m_rhoend;
        else
          m_rho *= Scalar( 0.1 );

        m_delta = max( m_delta, m_rho );  // Mantieni delta almeno rho

        // Output diagnostico
        if ( m_print_level >= 2 )
        {
          fmt::print(
              "\n"
              "    New RHO                   = {}\n"
              "    Number of function values = {}\n"
              "    Least value of F          = {}\n"
              "    The corresponding X is: {}\n",
              m_rho, m_num_f_eval, m_f_val( m_kopt ), print_vec( m_x_base + m_x_opt, 6 ) );
        }

        // Reset parametri per nuova fase
        ntrits        = 0;
        n_num_f_saved = m_num_f_eval;
        return Phase::TRUST_REGION;  // Continua con nuovo rho
      }
      if ( ntrits == -1 ) return Phase::EVALUATE;  // Caso speciale: nessuna valutazione
      return Phase::DONE;  // Convergenza raggiunta
    };

    Phase phase      = Phase::UPDATE_GRADIENT;
    bool  keep_going = true;

    if ( m_num_f_eval > m_max_f_eval )
    {
      m_reason = "CALFUN has been called MAXFUN times";
      m_status = Status::BOBYQA_TOO_MANY_EVALUATIONS;
      // finalization at bottom
      goto FINALIZE;
    }

    while ( keep_going )
    {
      switch ( phase )
      {
        case Phase::UPDATE_GRADIENT:
          info("UPDATE_GRADIENT");
          phase = phase_update_gradient();
          // std::cout << "m_knew = " << m_knew << '\n';
          break;

        case Phase::TRUST_REGION:
          info("TRUST_REGION");
          phase = phase_trust_region();
          // std::cout << "m_knew = " << m_knew << '\n';
          break;

        case Phase::SHIFT_BASE:
          info("SHIFT_BASE");
          phase = phase_shift_base();
          // std::cout << "m_knew = " << m_knew << '\n';
          break;

        case Phase::RESCUE:
          std::cout << "\\n";
          info("RESCUE");
          phase = phase_rescue();
          // std::cout << "m_knew = " << m_knew << '\n';
          break;

        case Phase::ALTMOV:
          info("ALTMOV");
          phase = phase_altmov();
          // std::cout << "m_knew = " << m_knew << '\n';
          break;

        case Phase::COMPUTE_VLAG:
          info("COMPUTE_VLAG");
          phase = phase_compute_vlag();
          // std::cout << "m_knew = " << m_knew << '\n';
          break;

        case Phase::EVALUATE:
          info("EVALUATE");
          phase = phase_evaluate();
          // std::cout << "m_knew = " << m_knew << '\n';
          break;

        case Phase::FIND_FAR:
          info("FIND_FAR");
          phase = phase_find_far();
          // std::cout << "m_knew = " << m_knew << '\n';
          break;

        case Phase::REDUCE_RHO:
          info("REDUCE_RHO");
          phase = phase_reduce_rho();
          // std::cout << "m_knew = " << m_knew << '\n';
          break;

        case Phase::DONE:
          info("DONE");
          //  Termina l'esecuzione dell'algoritmo dopo aver raggiunto convergenza
          //  o dopo aver esaurito le risorse di calcolo.
          keep_going = false;
          break;

        case Phase::ERROR:
          // std::cout << "ERROR m_knew = " << m_knew << '\n';
          //  Gestisce gli errori dell'algoritmo, impostando il messaggio di errore
          //  e lo stato appropriato.
          if ( m_print_level > 0 && m_reason.length() > 0 ) print_error( m_reason );
          keep_going = false;
          break;

      }  // switch
    }  // while

  FINALIZE:
    // finalize result: set X to best point if smaller than starting fsave
    if ( m_f_val( m_kopt ) <= fsave )
    {
      X = ( m_x_base + m_x_opt ).cwiseMax( m_x_lower ).cwiseMin( m_x_upper );
      f = m_f_val( m_kopt );
    }

    if ( m_print_level >= 1 )
    {
      fmt::print(
          "\n"
          "    At the return from BOBYQA\n"
          "    Number of function values = {}\n"
          "    Least value of F          = {}\n"
          "    The corresponding X is: {}\n",
          m_num_f_eval, f, print_vec( X, 6 ) );
    }

    // non ha senso lo rimuovo
    // if ( status == Status::BOBYQA_SUCCESS ) m_xbase( 0 ) = f;
    return m_status;
  }

  /*----------------------------------------------------------------------------*/
  template <typename Scalar>
  void
  BOBYQA_minimizer<Scalar>::prelim( Vector & X )
  {
    /* Set XBASE to the initial vector of variables, and set the initial elements
       of XPT, BMAT, HQ, PQ and ZMAT to zero. */
    m_xpt.setZero();
    m_B.setZero();
    m_HQ.setZero();
    m_pq.setZero();
    m_Z.setZero();

    m_x_base.noalias() = X.cwiseMin( m_x_upper ).cwiseMax( m_x_lower );

    auto OBJ = [&]( integer const ipos ) -> void
    {
      Scalar & f = m_f_val( ipos );
      f          = eval( X );
      if ( f < m_f_val( m_kopt ) ) m_kopt = ipos;
    };

    integer const & n  = m_nv;
    auto const &    su = m_s_upper;
    auto const &    sl = m_s_lower;

    Scalar  rhosq = m_rhobeg * m_rhobeg;
    Scalar  recip = Scalar(1) / rhosq;
    integer np    = n + 1;
    Scalar  fbeg  = 0;

    m_x_base = X;
    m_xpt.setZero();
    m_B.setZero();
    m_HQ.setZero();
    m_pq.setZero();
    m_Z.setZero();

    // =============== PARTE 1: PUNTO BASE (nf = 1) ===============
    {
      m_f_val( 0 ) = fbeg = eval( m_x_base );
      m_kopt              = 0;
    }

    // =============== PARTE 2: PUNTI 2..n+1 (perturbazioni positive) ===============
    for ( integer nf = 2; nf <= n + 1; ++nf )
    {
      integer point_idx = nf - 1;  // Indice del punto corrente in xpt
      integer coord_idx = nf - 2;  // Indice della coordinata su cui operiamo

      Scalar step = m_rhobeg;
      if ( abs( su( coord_idx ) ) < m_eps ) step = -step;

      // adatta passo nel caso vicino al bordo
      auto const & SL = sl( coord_idx );
      auto const & SU = su( coord_idx );
      if ( abs( step - SL ) < m_eps ) step = SL;
      if ( abs( step - SU ) < m_eps ) step = SU;

      m_xpt( coord_idx, point_idx ) = step;

      auto &   XX = X( coord_idx );
      Scalar & f  = m_f_val( point_idx );
      XX += step;
      OBJ( point_idx );
      XX -= step;

      // Aggiorna il miglior punto
      if ( f < m_f_val( m_kopt ) ) m_kopt = point_idx;

      // Aggiorna il modello quadratico per i punti 2..n+1
      m_g_opt( coord_idx ) = ( f - fbeg ) / step;

      if ( m_npt < nf + n )
      {
        m_B( coord_idx, 0 )                 = -Scalar(1) / step;
        m_B( coord_idx, point_idx )         = Scalar(1) / step;
        m_B( coord_idx, m_npt + point_idx ) = -Scalar(0.5) * rhosq;
      }
    }

    // =============== PARTE 3: PUNTI n+2..2n+1 (perturbazioni negative) ===============
    integer nf_max = min( 2 * n + 1, m_npt );
    for ( integer nf = n + 2; nf <= nf_max; ++nf )
    {
      integer point_idx       = nf - 1;      // Indice del punto corrente in xpt
      integer coord_idx       = nf - n - 2;  // Indice della coordinata su cui operiamo
      integer other_point_idx = nf - n - 1;  // Indice del punto corrispondente con perturbazione positiva

      Scalar step_pos = m_xpt( coord_idx, other_point_idx );  // Step positivo precedente
      Scalar step_neg = -m_rhobeg;                            // Step negativo di default

      auto const & SU = su( coord_idx );
      auto const & SL = sl( coord_idx );

      // Adatta stepb ai limiti
      if ( abs( SL ) < m_eps ) step_neg = min( Scalar(2) * m_rhobeg, SU );
      if ( abs( SU ) < m_eps ) step_neg = max( -Scalar(2) * m_rhobeg, SL );

      if ( abs( step_neg - SL ) < m_eps ) step_neg = SL;
      if ( abs( step_neg - SU ) < m_eps ) step_neg = SU;

      m_xpt( coord_idx, point_idx ) = step_neg;

      auto &   XX = X( coord_idx );
      Scalar & f  = m_f_val( point_idx );

      XX += step_neg;
      OBJ( point_idx );
      XX -= step_neg;

      // Aggiorna il modello quadratico
      Scalar temp = ( f - fbeg ) / step_neg;
      Scalar diff = step_neg - step_pos;

      auto & g                     = m_g_opt( coord_idx );
      m_HQ( coord_idx, coord_idx ) = Scalar(2) * ( temp - g ) / diff;
      g                            = ( g * step_neg - temp * step_pos ) / diff;

      // Scambia punti se migliora il modello
      Scalar & fo = m_f_val( other_point_idx );
      if ( step_pos * step_neg < 0 && f < fo )
      {
        std::swap( f, fo );
        std::swap( step_neg, step_pos );
        if ( m_kopt == point_idx ) m_kopt = other_point_idx;
        m_xpt( coord_idx, other_point_idx ) = step_pos;
        m_xpt( coord_idx, point_idx )       = step_neg;
      }

      // Aggiorna B e Z matrices
      Scalar t1                         = ( step_pos + step_neg ) / ( step_pos * step_neg );
      Scalar t2                         = Scalar(0.5) / step_pos;
      m_B( coord_idx, 0 )               = -t1;
      m_B( coord_idx, point_idx )       = -t2;
      m_B( coord_idx, other_point_idx ) = t1 + t2;

      Scalar t3                         = sqrt( Scalar(2) ) / ( step_pos * step_neg );
      Scalar t4                         = sqrt( Scalar(0.5) ) / rhosq;
      m_Z( 0, coord_idx )               = t3;
      m_Z( point_idx, coord_idx )       = t4;
      m_Z( other_point_idx, coord_idx ) = -( t3 + t4 );
    }

    // =============== PARTE 4: PUNTI 2n+2..max_iterations (combinazioni) ===============
    nf_max = min( 2 * n + 1, m_npt );
    for ( integer nf = 2 * np; nf <= nf_max; ++nf )
    {
      integer point_idx    = nf - 1;      // Indice del punto corrente in xpt
      integer zmat_col_idx = nf - n - 2;  // Indice della colonna in zmat

      // Calcola ipt e jpt usando la formula del BOBYQA originale
      integer temp_val = ( point_idx - np ) / n;
      integer jpt      = point_idx - temp_val * np - 1;
      integer ipt      = jpt + temp_val;

      if ( ipt >= n )
      {
        temp_val = jpt;
        jpt      = ipt - n;
        ipt      = temp_val;
      }

      // Crea nuovo punto combinando due punti esistenti
      // NOTA: usiamo ipt+1 e jpt+1 perché ipt e jpt sono indici di coordinate
      // e i punti corrispondenti hanno indici ipt+1 e jpt+1 (1-based in Fortran)
      integer point_i = ipt + 1;  // Punto con perturbazione positiva lungo ipt
      integer point_j = jpt + 1;  // Punto con perturbazione positiva lungo jpt

      // Crea nuovo punto combinando due punti esistenti
      Scalar stepi = m_xpt( ipt, point_i );
      Scalar stepj = m_xpt( jpt, point_j );

      // Valuta la funzione nel nuovo punto
      {
        auto const & SU = su( ipt );
        auto const & SL = sl( ipt );

        // Adatta stepb ai limiti
        if ( abs( SL ) < m_eps ) stepi = min( Scalar(2) * m_rhobeg, SU );
        if ( abs( SU ) < m_eps ) stepi = max( -Scalar(2) * m_rhobeg, SL );

        if ( abs( stepi - SL ) < m_eps ) stepi = SL;
        if ( abs( stepi - SU ) < m_eps ) stepi = SU;
      }
      {
        auto const & SU = su( jpt );
        auto const & SL = sl( jpt );

        // Adatta stepb ai limiti
        if ( abs( SL ) < m_eps ) stepj = min( Scalar(2) * m_rhobeg, SU );
        if ( abs( SU ) < m_eps ) stepj = max( -Scalar(2) * m_rhobeg, SL );

        if ( abs( stepj - SL ) < m_eps ) stepj = SL;
        if ( abs( stepj - SU ) < m_eps ) stepj = SU;
      }

      auto &   Xi = X( ipt );
      auto &   Xj = X( jpt );
      Scalar & f  = m_f_val( point_idx );
      Scalar & fi = m_f_val( point_i );
      Scalar & fj = m_f_val( point_j );

      Xi += stepi;
      Xj += stepj;
      OBJ( point_idx );
      Xi -= stepi;
      Xj -= stepj;

      // Crea nuovo punto combinando due punti esistenti
      m_xpt( ipt, point_idx ) = stepi;
      m_xpt( jpt, point_idx ) = stepj;

      // Aggiorna Z matrix e Hessiana
      // NOTA: anche in zmat usiamo ipt+1 e jpt+1
      m_Z( 0, zmat_col_idx )         = recip;
      m_Z( point_idx, zmat_col_idx ) = recip;
      m_Z( point_i, zmat_col_idx )   = -recip;
      m_Z( point_j, zmat_col_idx )   = -recip;

      Scalar product = stepi * stepj;
      // NOTA: anche qui usiamo fval[point_i] e fval[point_j]
      if ( ipt > jpt ) std::swap( ipt, jpt );
      m_HQ( ipt, jpt ) = ( fbeg + f - fi - fj ) / product;
    }

    // Aggiorna x al miglior punto trovato
    X = m_x_base + m_xpt.col( m_kopt );

  } /* prelim */


  template <typename Scalar>
  void
  BOBYQA_minimizer<Scalar>::altmov()
  {
    const Scalar one_plus_sqrt2 = 1 + sqrt( Scalar( 2 ) );

    //Vector W( 2 * m_nv );

    Scalar  csave  = 0;
    Scalar  stpsav = 0;
    Scalar  step   = 0;
    integer ksav   = -1;
    integer ibdsav = 0;

    // Compute H column for knew
    m_hcol.setZero();
    for ( integer j = 0; j < m_nptm; ++j ) m_hcol.noalias() += m_Z( m_knew, j ) * m_Z.col( j );
    m_alpha   = m_hcol( m_knew );
    Scalar ha = Scalar( 0.5 ) * m_alpha;

    // Compute gradient glag
    m_g_lag.noalias() = m_B.col( m_knew );

    for ( integer k = 0; k < m_npt; ++k )
    {
      auto const & xptk = m_xpt.col( k );
      Scalar       temp = xptk.dot( m_x_opt );
      temp *= m_hcol( k );
      m_g_lag.head( m_nv ) += temp * xptk;
    }

    // Search best line
    Scalar presav = 0;
    for ( integer k = 0; k < m_npt; ++k )
    {
      if ( k == m_kopt ) continue;

      Vector tmp = m_xpt.col( k ) - m_x_opt;
      Scalar dderiv = tmp.dot( m_g_lag.head( m_nv ) );
      m_distsq      = tmp.squaredNorm();

      Scalar  subd = m_adelt / sqrt( m_distsq );
      Scalar  slbd = -subd;
      integer ilbd = 0;
      integer iubd = 0;

      Scalar sumin = min( Scalar( 1 ), subd );

      // bound projection - CORRETTA come originale
      for ( integer i = 0; i < m_nv; ++i )
      {
        auto const & xo = m_x_opt( i );
        auto const & su = m_s_upper( i );
        auto const & sl = m_s_lower( i );

        Scalar temp = m_xpt( i, k ) - xo;

        if ( temp > 0 )
        {
          if ( slbd * temp < sl - xo )
          {
            slbd = ( sl - xo ) / temp;
            ilbd = -i - 1;
          }
          if ( subd * temp > su - xo )
          {
            subd = max( ( su - xo ) / temp, sumin );
            iubd = i + 1;
          }
        }
        else if ( temp < 0 )
        {
          if ( slbd * temp > su - xo )
          {
            slbd = ( su - xo ) / temp;
            ilbd = i + 1;
          }
          if ( subd * temp < sl - xo )
          {
            subd = max( ( sl - xo ) / temp, sumin );
            iubd = -i - 1;
          }
        }
      }

      Scalar  vlag;
      integer isbd;

      if ( k == m_knew )
      {
        Scalar diff = dderiv - 1;
        step        = slbd;
        vlag        = slbd * ( dderiv - slbd * diff );
        isbd        = ilbd;

        Scalar temp = subd * ( dderiv - subd * diff );
        if ( abs( temp ) > abs( vlag ) )
        {
          step = subd;
          vlag = temp;
          isbd = iubd;
        }

        Scalar tempd = Scalar( 0.5 ) * dderiv;
        Scalar tempa = tempd - diff * slbd;
        Scalar tempb = tempd - diff * subd;
        if ( tempa * tempb < 0 )
        {
          temp = tempd * tempd / diff;
          if ( abs( temp ) > abs( vlag ) )
          {
            step = tempd / diff;
            vlag = temp;
            isbd = 0;
          }
        }
      }
      else
      {
        step = slbd;
        vlag = slbd * ( 1 - slbd );
        isbd = ilbd;

        Scalar temp = subd * ( 1 - subd );
        if ( abs( temp ) > abs( vlag ) )
        {
          step = subd;
          vlag = temp;
          isbd = iubd;
        }

        if ( subd > Scalar( 0.5 ) )
        {
          if ( abs( vlag ) < Scalar( 0.25 ) )
          {
            step = Scalar( 0.5 );
            vlag = Scalar( 0.25 );
            isbd = 0;
          }
        }

        vlag *= dderiv;
      }

      Scalar temp   = step * ( 1 - step ) * m_distsq;
      Scalar predsq = power2( vlag ) * ( power2( vlag ) + ha * power2( temp ) );

      if ( predsq > presav )
      {
        presav = predsq;
        ksav   = k;
        stpsav = step;
        ibdsav = isbd;
      }
    }

    // construct xnew - Se ksav è ancora 0, c'è un problema
    if ( ksav == -1 )
    {
      // Fallback: usa knew come ksav
      ksav = m_knew;
      // Calcola un step di default
      Scalar dist = ( m_xpt.col( ksav ) - m_x_opt ).norm();
      stpsav      = m_adelt / dist;
      ibdsav      = 0;
    }

    m_x_new = ( m_x_opt + stpsav * ( m_xpt.col( ksav ) - m_x_opt ) ).cwiseMin( m_s_upper ).cwiseMax( m_s_lower );

    // Applica i bound specifici
    if ( ibdsav < 0 )
    {
      integer idx    = -ibdsav - 1;
      m_x_new( idx ) = m_s_lower( idx );
    }
    if ( ibdsav > 0 )
    {
      integer idx    = ibdsav - 1;
      m_x_new( idx ) = m_s_upper( idx );
    }

    // ====== Cauchy step evaluation ======
    auto compute_cauchy_step = [&]( bool flip_grad )
    {
      if ( flip_grad ) m_g_lag.noalias() = -m_g_lag;

      Scalar bigstp = m_adelt + m_adelt;
      Scalar wfixsq = 0;
      Scalar ggfree = 0;

      Vector W( m_nv );

      for ( integer i = 0; i < m_nv; ++i )
      {
        Scalar tempa = std::min( m_x_opt( i ) - m_s_lower( i ), m_g_lag( i ) );
        Scalar tempb = std::max( m_x_opt( i ) - m_s_upper( i ), m_g_lag( i ) );

        if ( tempa > 0 || tempb < 0 )
        {
          W( i ) = bigstp;
          ggfree += m_g_lag( i ) * m_g_lag( i );
        }
      }

      if ( is_zero( ggfree ) )
      {
        m_x_alt.setZero();
        return Scalar( 0 );
      }

      Scalar step_local = 0;
      // recheck bounds until stable
      while ( true )
      {
        Scalar temp = power2( m_adelt ) - wfixsq;
        if ( temp <= 0 ) break;

        Scalar wsqsav = wfixsq;
        step_local    = sqrt( temp / ggfree );
        ggfree        = 0;

        for ( integer i = 0; i < m_nv; ++i )
        {
          if ( W( i ) == bigstp )
          {
            Scalar cand = m_x_opt( i ) - step_local * m_g_lag( i );
            if ( cand <= m_s_lower( i ) )
            {
              W( i ) = m_s_lower( i ) - m_x_opt( i );
              wfixsq += W( i ) * W( i );
            }
            else if ( cand >= m_s_upper( i ) )
            {
              W( i ) = m_s_upper( i ) - m_x_opt( i );
              wfixsq += W( i ) * W( i );
            }
            else
            {
              ggfree += m_g_lag( i ) * m_g_lag( i );
            }
          }
        }
        if ( !( wfixsq > wsqsav && ggfree > 0 ) ) break;
      }

      Scalar gw = 0;
      for ( integer i = 0; i < m_nv; ++i )
      {
        if ( W( i ) == bigstp )
        {
          W( i )       = -step_local * m_g_lag( i );
          Scalar v     = std::min( m_x_opt( i ) + W( i ), m_s_upper( i ) );
          m_x_alt( i ) = std::max( v, m_s_lower( i ) );
        }
        else if ( is_zero( W( i ) ) )
          m_x_alt( i ) = m_x_opt( i );
        else if ( m_g_lag( i ) > 0 )
          m_x_alt( i ) = m_s_lower( i );
        else
          m_x_alt( i ) = m_s_upper( i );

        gw += m_g_lag( i ) * W( i );
      }

      Scalar curv = 0;
      for ( integer k = 0; k < m_npt; ++k )
      {
        Scalar temp = W.dot( m_xpt.col( k ) );
        curv += m_hcol( k ) * temp * temp;
      }

      if ( flip_grad ) curv = -curv;

      Scalar result;
      if ( curv > -gw && curv < -one_plus_sqrt2 * gw )
      {
        Scalar scale = -gw / curv;
        for ( integer i = 0; i < m_nv; ++i )
        {
          Scalar v     = std::min( m_x_opt( i ) + scale * W( i ), m_s_upper( i ) );
          m_x_alt( i ) = std::max( v, m_s_lower( i ) );
        }
        Scalar temp = Scalar( 0.5 ) * gw * scale;
        result      = temp * temp;
      }
      else
      {
        Scalar temp = gw + Scalar( 0.5 ) * curv;
        result      = temp * temp;
      }

      return result;
    };

    // Evaluate downhill and uphill version
    Scalar c1 = compute_cauchy_step( false );
    Vector X_ALT_SAVE = m_x_alt;
    csave = c1;

    Scalar c2 = compute_cauchy_step( true );

    if ( csave > c2 )
    {
      m_x_alt = X_ALT_SAVE;
      m_cauchy = csave;
    }
    else
    {
      m_cauchy = c2;
    }
  }

  /**
   * @brief Procedura di salvataggio (Rescue) per ripristinare la geometria dell'interpolazione.
   *
   * @details Questa funzione viene chiamata quando l'algoritmo rileva instabilità numerica
   * o quando i punti di interpolazione sono degenerati (es. quasi allineati).
   * La procedura:
   * 1. Sposta l'origine delle coordinate nel punto ottimo corrente (XOPT).
   * 2. Ripristina la matrice Z (che gestisce l'indipendenza lineare).
   * 3. Aggiorna la matrice Hessiana (HQ) e le matrici di lavoro (BMAT).
   * 4. Tenta di sostituire iterativamente i punti peggiori con punti geometricamente
   * più validi (basati sui bounds o sulla direzione di massima curvatura),
   * utilizzando operazioni vettoriali per efficienza.
   */


  template <typename Scalar>
  void
  BOBYQA_minimizer<Scalar>::rescue()
  {
    // Vettore di lavoro locale: prime m_dim entrate per calcoli su variabili,
    // successive m_npt per calcoli sui punti.
    Vector m_work( m_dim + m_npt );

    // Salva FVAL(KOPT) prima dello shift (come FBASE in Fortran)
    Scalar fbase = m_f_val( m_kopt );

    m_beta       = 0;
    Scalar sumpq = 0;
    Scalar winc  = 0;

    const integer np    = m_nv + 1;  // N+1
    const Scalar  sfrac = Scalar( 0.5 ) / Scalar( np );

    // ==============================================================================
    // STEP 1: Shift delle coordinate all'origine e calcolo distanze
    // ==============================================================================
    // Spostiamo tutti i punti XPT in modo che XOPT diventi l'origine (0,0...).
    // XPT è (n x npt) in C++ vs (npt x n) in Fortran
    for ( integer k = 0; k < m_npt; ++k ) { m_xpt.col( k ) -= m_x_opt; }

    // Calcoliamo la somma dei pesi PQ
    sumpq = m_pq.sum();

    // Calcoliamo le norme quadrate di ogni punto dopo lo shift
    for ( integer k = 0; k < m_npt; ++k ) { m_work( m_dim + k ) = m_xpt.col( k ).squaredNorm(); }

    // Troviamo la distanza massima per il calcolo di winc
    winc = m_work.tail( m_npt ).maxCoeff();

    // Resettiamo ZMAT: le prime m_nptm colonne vengono azzerate
    if ( m_nptm > 0 ) { m_Z.leftCols( m_nptm ).setZero(); }

    // ==============================================================================
    // STEP 2: Aggiornamento Hessiana (HQ) dovuto allo shift
    // ==============================================================================
    // Calcoliamo il vettore di correzione: v = 0.5 * sumpq * XOPT + XPT^T * PQ
    Vector vec_correction = ( Scalar( 0.5 ) * sumpq ) * m_x_opt;
    // Nota: m_xpt è (n x npt), m_pq è (npt), quindi m_xpt * m_pq = vettore n x 1
    vec_correction += m_xpt * m_pq;

    // Aggiornamento packed upper-triangular di HQ
    for ( integer j = 0; j < m_nv; ++j )
    {
      Scalar wj     = vec_correction( j );
      Scalar xopt_j = m_x_opt( j );
      for ( integer i = 0; i <= j; ++i ) { m_HQ( i, j ) += vec_correction( i ) * xopt_j + wj * m_x_opt( i ); }
    }

    // ==============================================================================
    // STEP 3: Shift variabili base e costruzione PTSAUX
    // ==============================================================================
    m_x_base += m_x_opt;
    m_s_lower -= m_x_opt;
    m_s_upper -= m_x_opt;

    // Costruiamo i punti ausiliari PTSAUX
    for ( integer j = 0; j < m_nv; ++j )
    {
      Scalar sl_j = m_s_lower( j );
      Scalar su_j = m_s_upper( j );

      Scalar cand0 = std::min( m_delta, su_j );
      Scalar cand1 = std::max( -m_delta, sl_j );

      // Ordinamento per garantire stabilità numerica
      if ( cand0 + cand1 < 0 ) std::swap( cand0, cand1 );

      // Riduzione se troppo vicino a zero
      if ( std::abs( cand1 ) < Scalar( 0.5 ) * std::abs( cand0 ) ) { cand1 = Scalar( 0.5 ) * cand0; }

      m_ptsaux( 0, j ) = cand0;
      m_ptsaux( 1, j ) = cand1;
    }

    // XOPT ora è formalmente l'origine
    m_x_opt.setZero();

    // Azzeriamo BMAT (Fortran: BMAT(I,J)=ZERO per tutte le I,J). In C++ m_B
    // è (N x NDIM) e pertanto dobbiamo azzerare l'intera matrice per avere
    // lo stesso effetto dell'implementazione Fortran.
    m_B.setZero();

    // ==============================================================================
    // STEP 4: Costruzione iniziale PTSAUX / BMAT / ZMAT
    // ==============================================================================
    m_ptsid( 0 ) = sfrac;  // Punto base (indice 0)

    for ( integer j = 0; j < m_nv; ++j )
    {
      // Indici 0-based per C++
      integer jp  = j + 1;      // Indice punto positivo (1-based in Fortran, 0-based in C++)
      integer jpn = jp + m_nv;  // Indice punto negativo (1-based in Fortran, 0-based in C++)

      m_ptsid( jp ) = Scalar( j + 1 ) + sfrac;  // j+1 perché in Fortran j parte da 1

      Scalar d0 = m_ptsaux( 0, j );
      Scalar d1 = m_ptsaux( 1, j );

      Scalar inv_d0   = Scalar( 1 ) / d0;
      Scalar diff_inv = Scalar( 1 ) / ( d0 - d1 );

      if ( jpn < m_npt )
      {
        m_ptsid( jpn ) = Scalar( j + 1 ) / Scalar( np ) + sfrac;

        // CORREZIONE: m_bmat è (n, m_dim) quindi:
        // - Prima coordinata: indice variabile (j)
        // - Seconda coordinata: indice punto (jp, jpn, 0)
        m_B( j, jp )  = -diff_inv + inv_d0;
        m_B( j, jpn ) = diff_inv + Scalar( 1 ) / d1;
        m_B( j, 0 )   = -m_B( j, jp ) - m_B( j, jpn );

        // ZMAT: (npt, m_nptm) come in Fortran
        Scalar z_base = std::sqrt( Scalar( 2.0 ) ) / std::abs( d0 * d1 );
        m_Z( 0, j )   = z_base;
        m_Z( jp, j )  = z_base * d1 * diff_inv;
        m_Z( jpn, j ) = -z_base * d0 * diff_inv;
      }
      else
      {
        // Caso in cui abbiamo meno punti
        m_B( j, 0 )  = -inv_d0;
        m_B( j, jp ) = inv_d0;
        // Termine quadratico: colonna j + npt (ma verifica che esista)
        if ( j + m_npt < m_dim ) { m_B( j, j + m_npt ) = -Scalar( 0.5 ) * power2( d0 ); }
      }
    }

    // Gestione punti extra (oltre 2*n + 1)
    if ( m_npt >= m_nv + np )
    {
      for ( integer k = 2 * np; k <= m_npt; ++k )
      {
        integer k_idx = k - 1;  // 0-based
        // Calcoli come in Fortran
        integer iw = static_cast<integer>( ( Scalar( k - np ) - Scalar( 0.5 ) ) / Scalar( m_nv ) );
        integer ip = k - np - iw * m_nv;
        integer iq = ip + iw;
        if ( iq > m_nv ) iq -= m_nv;

        m_ptsid( k_idx ) = Scalar( ip ) + Scalar( iq ) / Scalar( np ) + sfrac;

        // Indici 0-based per array
        integer ip0 = ip - 1;
        integer iq0 = iq - 1;

        Scalar tmp = Scalar( 1 ) / ( m_ptsaux( 0, ip0 ) * m_ptsaux( 0, iq0 ) );

        integer col_z       = k - np - 1;
        m_Z( 0, col_z )     = tmp;
        m_Z( ip0, col_z )   = -tmp;
        m_Z( iq0, col_z )   = -tmp;
        m_Z( k_idx, col_z ) = tmp;
      }
    }

    // ==============================================================================
    // MAIN PHASE LOOP: Ripristino iterativo
    // ==============================================================================
    integer nrem = m_npt;
    integer kold = 0;
    m_knew       = m_kopt;

    // Funzione helper per lo swap - CORRETTA per indici 0-based
    auto swap_points = [&]( integer k1, integer k2 )
    {
      if ( k1 == k2 ) return;

      // Scambia colonne di BMAT (perché m_bmat è (n, m_dim) e i punti sono colonne)
      m_B.col( k1 ).swap( m_B.col( k2 ) );

      // Scambia righe di ZMAT (perché m_zmat è (npt, m_nptm) e i punti sono righe)
      if ( m_nptm > 0 ) { m_Z.row( k1 ).swap( m_Z.row( k2 ) ); }

      // Scambia altri array
      std::swap( m_ptsid( k1 ), m_ptsid( k2 ) );
      std::swap( m_work( m_dim + k1 ), m_work( m_dim + k2 ) );
      std::swap( m_v_lag( k1 ), m_v_lag( k2 ) );
    };

    while ( nrem > 0 )
    {
      // PHASE 1: Scambia kold e knew
      if ( m_knew != m_kopt )
      {
        swap_points( kold, m_knew );
        m_work( m_dim + m_knew ) = 0;
        nrem--;

        // Aggiorna BMAT e ZMAT
        update();  // Assicurati che update() sia corretta per indici 0-based

        if ( nrem == 0 ) break;

        m_work.tail( m_npt ) = m_work.tail( m_npt ).cwiseAbs();
      }

      // PHASE 2: Trova nuovo candidato knew
      Scalar  dsqmin    = 0;
      integer best_knew = -1;

      for ( integer k = 0; k < m_npt; ++k )
      {
        Scalar d = m_work( m_dim + k );
        if ( d > 0 && ( is_zero( dsqmin ) || d < dsqmin ) )
        {
          dsqmin    = d;
          best_knew = k;
        }
      }

      if ( best_knew == -1 ) break;
      m_knew = best_knew;

      // PHASE 3: Calcola w_vec (come nel codice originale)
      Vector w_vec( m_npt + m_nv );
      w_vec.tail( m_nv ) = m_xpt.col( m_knew );

      for ( integer k = 0; k < m_npt; ++k )
      {
        if ( k == m_kopt )
        {
          w_vec( k ) = 0;
          continue;
        }

        Scalar sum = 0;
        if ( is_zero( m_ptsid( k ) ) )
        {
          // Punto originale: prodotto scalare
          sum = m_xpt.col( k ).dot( m_xpt.col( m_knew ) );
        }
        else
        {
          // Punto artificiale - DECODIFICA CORRETTA
          // In Fortran: IP = PTSID(K), IQ = (PTSID(K)-IP) * (N+1)

          integer ip      = static_cast<integer>( m_ptsid( k ) );
          Scalar  iq_real = np * m_ptsid( k ) - ip * np;
          integer iq      = static_cast<integer>( iq_real );  // troncamento

          if ( ip > 0 ) { sum += w_vec( m_npt + ip - 1 ) * m_ptsaux( 0, ip - 1 ); }
          if ( iq > 0 )
          {
            integer iw = 1;         // Usa PTSAUX(1, iq) per default
            if ( ip == 0 ) iw = 2;  // Usa PTSAUX(2, iq) se ip = 0
            sum += w_vec( m_npt + iq - 1 ) * m_ptsaux( iw - 1, iq - 1 );
          }
        }
        w_vec( k ) = Scalar( 0.5 ) * power2( sum );
      }

      // PHASE 4: Calcola VLAG e BETA
      // Parte BMAT: VLAG(K) = Σ_J BMAT(K,J) * W(NPT+J)
      // In C++: m_bmat è (n, m_dim), quindi m_bmat.col(k) è il punto k
      for ( integer k = 0; k < m_npt; ++k ) { m_v_lag( k ) = m_B.col( k ).dot( w_vec.tail( m_nv ) ); }

      // Parte ZMAT
      m_beta = 0;
      for ( integer j = 0; j < m_nptm; ++j )
      {
        Scalar sum_z = m_Z.col( j ).dot( w_vec.head( m_npt ) );
        m_beta -= power2( sum_z );

        // Aggiorna VLAG
        m_v_lag.head( m_npt ) += sum_z * m_Z.col( j );
      }

      // Calcola bsum e distsq
      Scalar bsum = 0;
      m_distsq    = m_xpt.col( m_knew ).squaredNorm();

      for ( integer j = 0; j < m_nv; ++j )
      {
        // Somma sulle prime npt colonne di BMAT (punti)
        Scalar sum = m_B.row( j ).dot( w_vec );
        bsum += sum * w_vec( m_npt + j );  // UNA SOLA VOLTA!
        m_v_lag( m_npt + j ) = sum;
      }

      m_beta = Scalar( 0.5 ) * power2( m_distsq ) + m_beta - bsum;
      m_v_lag( m_kopt ) += Scalar( 1.0 );

      // PHASE 5: Trova kold con il massimo denominatore
      Scalar vlmxsq = 0;
      for ( integer k = 0; k < m_npt; ++k )
      {
        Scalar temp = power2( m_v_lag( k ) );
        if ( temp > vlmxsq ) vlmxsq = temp;
      }

      Scalar  denom_max = 0;
      integer best_kold = -1;

      for ( integer k = 0; k < m_npt; ++k )
      {
        if ( !is_zero( m_ptsid( k ) ) )
        {
          Scalar hdiag = 0;
          if ( m_nptm > 0 )
          {
            auto R = m_Z.row( k );
            hdiag  = R.dot( R );
          }
          Scalar den = m_beta * hdiag + power2( m_v_lag( k ) );
          if ( den > denom_max )
          {
            denom_max = den;
            best_kold = k;
          }
        }
      }

      if ( best_kold == -1 ) break;
      kold = best_kold;

      // PHASE 6: Verifica denominatore
      if ( denom_max <= Scalar( 0.01 ) * vlmxsq )
      {
        m_work( m_dim + m_knew ) = -m_work( m_dim + m_knew ) - winc;
        continue;
      }

      // PHASE 7: Gestione caso speciale (knew == kopt)
      if ( m_knew == m_kopt )
      {
        m_work( m_dim + m_knew ) = 0;
        nrem--;
        if ( nrem == 0 ) break;
      }
    }

    // ==============================================================================
    // FASE FINALE: Valutazione e aggiornamento modello (L260-L350 in Fortran)
    // ==============================================================================
    for ( integer kpt = 0; kpt < m_npt; ++kpt )
    {
      if ( is_zero( m_ptsid( kpt ) ) ) continue;

      // Controllo limite valutazioni funzione
      if ( m_num_f_eval >= m_max_f_eval ) return;

      // STEP 1: Calcolo VQUAD (valore del modello quadratico corrente)
      Scalar vquad = fbase;

      // Decodifica PTSID
      integer ip      = static_cast<integer>( m_ptsid( kpt ) );
      Scalar  iq_real = np * m_ptsid( kpt ) - ip * np;
      integer iq      = static_cast<integer>( iq_real );  // troncamento

      Scalar xp = 0, xq = 0;

      if ( ip > 0 && ip <= m_nv )
      {
        xp = m_ptsaux( 0, ip - 1 );  // -1 per 0-based
        vquad += xp * ( m_g_opt( ip - 1 ) + Scalar( 0.5 ) * xp * m_HQ( ip - 1, ip - 1 ) );
      }

      if ( iq > 0 && iq <= m_nv )
      {
        if ( ip == 0 ) { xq = m_ptsaux( 1, iq - 1 ); }
        else
        {
          xq = m_ptsaux( 0, iq - 1 );
        }
        vquad += xq * ( m_g_opt( iq - 1 ) + Scalar( 0.5 ) * xq * m_HQ( iq - 1, iq - 1 ) );

        if ( ip > 0 && iq > 0 )
        {
          if ( ip <= iq )
            vquad += xp * xq * m_HQ( ip - 1, iq - 1 );
          else
            vquad += xp * xq * m_HQ( iq - 1, ip - 1 );
        }
      }

      // Write back new coordinates into m_xpt as the Fortran code does
      // (XPT(KPT,IP)=XP and XPT(KPT,IQ)=XQ) so that subsequent model
      // evaluation and CALFUN use the correct coordinates.
      if ( ip > 0 && ip <= m_nv ) { m_xpt( ip - 1, kpt ) = xp; }
      if ( iq > 0 && iq <= m_nv ) { m_xpt( iq - 1, kpt ) = xq; }

      // Contributo termini quadratici PQ
      for ( integer k = 0; k < m_npt; ++k )
      {
        Scalar temp = 0;
        if ( ip > 0 && ip <= m_nv )
        {
          // m_xpt(ip-1, k) = componente (ip-1) del punto k
          temp += xp * m_xpt( ip - 1, k );
        }
        if ( iq > 0 && iq <= m_nv ) { temp += xq * m_xpt( iq - 1, k ); }
        vquad += Scalar( 0.5 ) * m_pq( k ) * temp * temp;
      }

      // STEP 2: Valutazione funzione con gestione bounds
      Vector w_eval( m_nv );
      w_eval = ( m_x_base + m_xpt.col( kpt ) ).cwiseMax( m_x_lower ).cwiseMin( m_x_upper );

      // Valuta funzione
      Scalar f_value = eval( w_eval );
      m_f_val( kpt ) = f_value;

      // Aggiorna KOPT se trovato minimo migliore
      if ( f_value < m_f_val( m_kopt ) ) { m_kopt = kpt; }

      // STEP 3: Aggiornamento modello quadratico
      Scalar diff = f_value - vquad;

      // Aggiorna GOPT: GOPT(I) += DIFF * BMAT(KPT, I)
      m_g_opt += diff * m_B.col( kpt );  // m_bmat(i, kpt) non m_bmat(kpt, i)!

      // Aggiorna HQ e PQ
      for ( integer k = 0; k < m_npt; ++k )
      {
        if ( is_zero( m_ptsid( k ) ) )
        {
          // Punto originale: aggiorna PQ
          Scalar sum_z = m_Z.row( k ).dot( m_Z.row( kpt ) );
          m_pq( k ) += diff * sum_z;
        }
        else
        {
          // Punto artificiale: aggiorna HQ
          integer ipk      = static_cast<integer>( m_ptsid( k ) );
          Scalar  iq_realk = np * m_ptsid( k ) - ipk * np;  // CORRETTO!
          integer iqk      = static_cast<integer>( iq_realk );

          Scalar sum_z = m_Z.row( k ).dot( m_Z.row( kpt ) );
          Scalar temp  = diff * sum_z;

          if ( ipk == 0 && iqk > 0 && iqk <= m_nv )
          {
            // Solo iq
            Scalar xq_val = m_ptsaux( 1, iqk - 1 );
            m_HQ( iqk - 1, iqk - 1 ) += temp * xq_val * xq_val;
          }
          else if ( ipk > 0 && ipk <= m_nv )
          {
            // ip e possibilmente iq
            Scalar xp_val = m_ptsaux( 0, ipk - 1 );
            m_HQ( ipk - 1, ipk - 1 ) += temp * xp_val * xp_val;

            if ( iqk > 0 && iqk <= m_nv )
            {
              Scalar xq_val = m_ptsaux( 0, iqk - 1 );
              m_HQ( iqk - 1, iqk - 1 ) += temp * xq_val * xq_val;

              // Termine incrociato
              if ( ipk <= iqk )
                m_HQ( ipk - 1, iqk - 1 ) += temp * xp_val * xq_val;
              else
                m_HQ( iqk - 1, ipk - 1 ) += temp * xp_val * xq_val;
            }
          }
        }
      }

      // Marca punto come processato
      m_ptsid( kpt ) = 0;
    }
  }

  /**
   * @file trsbox.hpp
   * @brief Trust region subproblem solver for BOBYQA algorithm
   * @details Implements the truncated conjugate gradient method with bound constraints
   * @author Translated from Powell's Fortran implementation
   * @date 2023
   */

  /**
   * @brief Solves the bound-constrained trust region subproblem for BOBYQA
   *
   * @details
   * This function implements a truncated conjugate gradient method with bound
   * constraints to solve the trust region subproblem:
   * \f[
   * \begin{aligned}
   * \min_{d} \quad & g^T d + \frac{1}{2} d^T H d \\
   * \text{s.t.} \quad & \|d\| \leq \Delta \\
   *                   & l \leq x_{\text{opt}} + d \leq u
   * \end{aligned}
   * \f]
   *
   * The algorithm consists of three main phases:
   * 1. Initialization and identification of active bounds
   * 2. Truncated conjugate gradient iterations with projections
   * 3. Alternative iterations on the trust region boundary
   *
   * @tparam Scalar Numeric type (float, double, long double)
   *
   * @note This is a direct translation of Powell's Fortran implementation
   *       with indices shifted from 1-based to 0-based
   * @note Matrices are transposed compared to Fortran: m_xpt is N×NPT (Fortran is NPT×N)
   * @see Powell, M.J.D., "The BOBYQA algorithm for bound constrained optimization
   *      without derivatives", 2009
   */
  template <typename Scalar>
  void
  BOBYQA_minimizer<Scalar>::trsbox()
  {
    // ============================================================================
    // DEFINIZIONE DEGLI STATI
    // ============================================================================
    enum class State
    {
      TCG_Init,          // Inizializzazione TCG (Blocco 20)
      TCG_Loop,          // Calcolo direzione s (Blocco 30)
      TCG_Steplength,    // Calcolo lunghezza passo (Blocco 50)
      TCG_Update,        // Aggiornamento e decisione transizione (Blocco 70-80)
      Boundary_Init,     // Inizializzazione Boundary Phase (Blocco 90)
      Boundary_Loop,     // Loop Boundary (Blocco 100)
      Alternating_Loop,  // Loop alternativo (Blocco 120)
      Final,             // Uscita (Blocco 190)
      Exit
    };

    State current_state = State::TCG_Init;

    // ============================================================================
    // CONTESTO DATI (Variabili condivise tra gli stati)
    // ============================================================================

    integer iterc = 0, nact = 0, itermax = 0, iact = 0, itcsav = 0;
    Scalar  delsq = 0, qred = 0;
    Scalar  beta = 0, ggsav = 0, gredsq = 0, sdec = 0;
    Scalar  blen = 0, ds = 0, resid = 0, shs = 0, stepsq = 0, stplen = 0;
    Scalar  dredsq = 0, dredg = 0, sredg = 0;
    Scalar  dhd = 0, dhs = 0;
    Scalar  angbd = 0, tempa = 0, tempb = 0, ratio = 0;
    integer isav = 0, iu = 0, xsav = 0;
    Scalar  redmax = 0, redsav = 0, rdprev = 0, rdnext = 0, rednew = 0;
    Scalar  angt = 0, cth = 0, sth = 0;

    // ============================================================================
    // INIZIALIZZAZIONE (Blocco 10)
    // ============================================================================
    iterc = 0;
    nact  = 0;

    m_xbdi.setZero();
    m_d.setZero();
    m_g_new = m_g_opt;

    // Setup bound indicators
    for ( integer i = 0; i < m_nv; ++i )
    {
      if ( m_x_opt( i ) <= m_s_lower( i ) && m_g_opt( i ) >= 0 )
      {
        m_xbdi( i ) = -1;
        ++nact;
      }
      else if ( m_x_opt( i ) >= m_s_upper( i ) && m_g_opt( i ) <= 0 )
      {
        m_xbdi( i ) = 1;
        ++nact;
      }
    }

    delsq    = m_delta * m_delta;
    qred     = 0;
    m_crvmin = Scalar(-1);

    // ============================================================================
    // MACCHINA A STATI PRINCIPALE
    // ============================================================================
    while ( current_state != State::Exit )
    {
      switch ( current_state )
      {
        // ======================================================================
        case State::TCG_Init:
          // ======================================================================
          {
            // Inizializzazione o restart del TCG (Blocco 20)
            beta          = 0;
            current_state = State::TCG_Loop;
            break;
          }

        // ======================================================================
        case State::TCG_Loop:
          // ======================================================================
          {
            // Calcola la direzione di ricerca S (Blocco 30)
            // s = beta * s - g_new, ma azzera dove xbdi != 0
            m_s = beta * m_s - m_g_new;

            // Azzera componenti fisse
            for ( integer i = 0; i < m_nv; ++i )
            {
              if ( m_xbdi( i ) != 0 ) m_s( i ) = 0;
            }

            stepsq = m_s.squaredNorm();

            if ( stepsq == 0 )
            {
              current_state = State::Final;
              break;
            }

            if ( beta == 0 )
            {
              gredsq  = stepsq;
              itermax = iterc + m_nv - nact;
            }

            if ( gredsq * delsq <= Scalar(1.0e-4) * qred * qred )
            {
              current_state = State::Final;
              break;
            }

            // Calcola hs = H * s usando matrice piena m_hq (upper triangular)
            m_hs.setZero();

            // Parte diagonale
            for ( integer i = 0; i < m_nv; ++i ) { m_hs( i ) += m_HQ( i, i ) * m_s( i ); }

            // Parte triangolare superiore
            for ( integer i = 0; i < m_nv - 1; ++i )
            {
              for ( integer j = i + 1; j < m_nv; ++j )
              {
                m_hs( i ) += m_HQ( i, j ) * m_s( j );
                m_hs( j ) += m_HQ( i, j ) * m_s( i );
              }
            }

            // Contributo da m_pq e m_xpt
            for ( integer k = 0; k < m_npt; ++k )
            {
              if ( m_pq( k ) != 0 )
              {
                Scalar temp_qts = m_pq( k ) * m_xpt.col( k ).dot( m_s );
                m_hs += temp_qts * m_xpt.col( k );
              }
            }

            current_state = State::TCG_Steplength;
            break;
          }

        // ======================================================================
        case State::TCG_Steplength:
          // ======================================================================
          {
            // Calcolo Steplength (Blocco 50)
            resid = delsq;
            ds    = 0;
            shs   = 0;

            for ( integer i = 0; i < m_nv; ++i )
            {
              if ( m_xbdi( i ) == 0 )
              {
                resid -= m_d( i ) * m_d( i );
                ds += m_s( i ) * m_d( i );
                shs += m_s( i ) * m_hs( i );
              }
            }

            if ( resid <= 0 )
            {
              current_state = State::Boundary_Init;
              break;
            }

            Scalar temp = sqrt( stepsq * resid + ds * ds );
            if ( ds < 0 ) { blen = ( temp - ds ) / stepsq; }
            else
            {
              blen = resid / ( temp + ds );
            }

            stplen = blen;
            if ( shs > 0 ) { stplen = std::min( blen, gredsq / shs ); }

            // Gestione dei limiti semplici (Blocco 60)
            iact = 0;
            for ( integer i = 0; i < m_nv; ++i )
            {
              if ( m_s( i ) != 0 )
              {
                Scalar xsum = m_x_opt( i ) + m_d( i );
                Scalar temp_bound;

                if ( m_s( i ) > 0 ) { temp_bound = ( m_s_upper( i ) - xsum ) / m_s( i ); }
                else
                {
                  temp_bound = ( m_s_lower( i ) - xsum ) / m_s( i );
                }

                if ( temp_bound < stplen )
                {
                  stplen = temp_bound;
                  iact   = i + 1;  // +1 per usare 0 come "nessuna variabile"
                }
              }
            }

            current_state = State::TCG_Update;
            break;
          }

        // ======================================================================
        case State::TCG_Update:
          // ======================================================================
          {
            // Aggiornamento di CRVMIN, GNEW, D e QRED (Blocco 70)
            sdec = 0;

            if ( stplen > 0 )
            {
              ++iterc;
              Scalar temp = shs / stepsq;

              if ( iact == 0 && temp > 0 )
              {
                if ( m_crvmin == Scalar(-1) )
                  m_crvmin = temp;
                else
                  m_crvmin = std::min( m_crvmin, temp );
              }

              ggsav  = gredsq;
              gredsq = 0;

              // Aggiorna g_new e d usando operazioni vettoriali
              m_g_new += stplen * m_hs;
              m_d += stplen * m_s;

              // Calcola gredsq solo per componenti libere
              for ( integer i = 0; i < m_nv; ++i )
              {
                if ( m_xbdi( i ) == 0 ) gredsq += m_g_new( i ) * m_g_new( i );
              }

              sdec = std::max( stplen * ( ggsav - Scalar( 0.5 ) * stplen * shs ), Scalar( 0 ) );
              qred += sdec;
            }

            // Decide la transizione (Blocco 80)
            if ( iact > 0 )
            {
              // Restart TCG: Fissa la variabile e ricomincia
              ++nact;
              integer idx   = iact - 1;
              m_xbdi( idx ) = m_s( idx ) < 0 ? -1: 1;
              delsq -= power2( m_d( idx ) );

              if ( delsq <= 0 ) { current_state = State::Boundary_Init; }
              else
              {
                current_state = State::TCG_Init;
              }
            }
            else if ( stplen < blen )
            {
              // Continuazione TCG
              if ( iterc == itermax ) { current_state = State::Final; }
              else if ( sdec <= 0.01 * qred ) { current_state = State::Final; }
              else
              {
                beta          = gredsq / ggsav;
                current_state = State::TCG_Loop;
              }
            }
            else
            {
              // TCG terminato: passa alla Boundary Phase
              current_state = State::Boundary_Init;
            }
            break;
          }

        // ======================================================================
        case State::Boundary_Init:
          // ======================================================================
          {
            // Inizializzazione Boundary Phase (Blocco 90)
            m_crvmin = 0;

            if ( nact >= m_nv - 1 )
            {
              current_state = State::Final;
              break;
            }

            current_state = State::Boundary_Loop;
            break;
          }

        // ======================================================================
        case State::Boundary_Loop:
          // ======================================================================
          {
            // Loop dell'iterazione Boundary (Blocco 100)
            if ( nact >= m_nv - 1 )
            {
              current_state = State::Final;
              break;
            }

            // Calcola scalari
            dredsq = 0;
            dredg  = 0;
            gredsq = 0;

            for ( integer i = 0; i < m_nv; ++i )
            {
              if ( m_xbdi( i ) == 0 )
              {
                dredsq += m_d( i ) * m_d( i );
                dredg += m_d( i ) * m_g_new( i );
                gredsq += m_g_new( i ) * m_g_new( i );
                m_s( i ) = m_d( i );
              }
              else
              {
                m_s( i ) = 0;
              }
            }

            itcsav = iterc;

            // Calcola hs = H * s
            m_hs.setZero();

            // Parte diagonale
            for ( integer i = 0; i < m_nv; ++i ) { m_hs( i ) += m_HQ( i, i ) * m_s( i ); }

            // Parte triangolare superiore
            for ( integer i = 0; i < m_nv - 1; ++i )
            {
              for ( integer j = i + 1; j < m_nv; ++j )
              {
                m_hs( i ) += m_HQ( i, j ) * m_s( j );
                m_hs( j ) += m_HQ( i, j ) * m_s( i );
              }
            }

            // Contributo da m_pq e m_xpt
            for ( integer k = 0; k < m_npt; ++k )
            {
              if ( m_pq( k ) != 0 )
              {
                Scalar temp_qts = m_pq( k ) * m_xpt.col( k ).dot( m_s );
                m_hs += temp_qts * m_xpt.col( k );
              }
            }

            // Inizializza HRED
            if ( iterc == itcsav ) { m_hred = m_hs; }

            current_state = State::Alternating_Loop;
            break;
          }

        // ======================================================================
        case State::Alternating_Loop:
          // ======================================================================
          {
            // Loop alternativo (Blocco 120)
            ++iterc;
            Scalar temp = gredsq * dredsq - dredg * dredg;

            if ( temp <= 1.0e-4 * qred * qred )
            {
              current_state = State::Final;
              break;
            }

            // Calcola la nuova direzione S
            temp = std::sqrt( temp );
            for ( integer i = 0; i < m_nv; ++i )
            {
              if ( m_xbdi( i ) == 0 ) { m_s( i ) = ( dredg * m_d( i ) - dredsq * m_g_new( i ) ) / temp; }
              else
              {
                m_s( i ) = 0;
              }
            }
            sredg = -temp;

            // Calcola ANGBD e controlla i limiti (Blocco 130)
            angbd             = 1;
            iact              = 0;
            bool boundary_hit = false;

            for ( integer i = 0; i < m_nv; ++i )
            {
              if ( m_xbdi( i ) == 0 )
              {
                tempa = m_x_opt( i ) + m_d( i ) - m_s_lower( i );
                tempb = m_s_upper( i ) - m_x_opt( i ) - m_d( i );

                if ( tempa <= 0 )
                {
                  ++nact;
                  m_xbdi( i )  = -1;
                  boundary_hit = true;
                  break;
                }
                else if ( tempb <= 0 )
                {
                  ++nact;
                  m_xbdi( i )  = 1;
                  boundary_hit = true;
                  break;
                }

                ratio      = 1;
                Scalar ssq = m_d( i ) * m_d( i ) + m_s( i ) * m_s( i );

                // Controlla il limite inferiore
                Scalar temp_bound = ssq - ( m_x_opt( i ) - m_s_lower( i ) ) * ( m_x_opt( i ) - m_s_lower( i ) );
                if ( temp_bound > 0 )
                {
                  temp_bound = std::sqrt( temp_bound ) - m_s( i );
                  if ( angbd * temp_bound > tempa )
                  {
                    angbd = tempa / temp_bound;
                    iact  = i + 1;
                    xsav  = -1;
                  }
                }

                // Controlla il limite superiore
                temp_bound = ssq - ( m_s_upper( i ) - m_x_opt( i ) ) * ( m_s_upper( i ) - m_x_opt( i ) );
                if ( temp_bound > 0 )
                {
                  temp_bound = std::sqrt( temp_bound ) + m_s( i );
                  if ( angbd * temp_bound > tempb )
                  {
                    angbd = tempb / temp_bound;
                    iact  = i + 1;
                    xsav  = 1;
                  }
                }
              }
            }

            if ( boundary_hit )
            {
              current_state = State::Boundary_Loop;
              break;
            }

            // Calcola hs = H * s
            m_hs.setZero();

            // Parte diagonale
            for ( integer i = 0; i < m_nv; ++i ) { m_hs( i ) += m_HQ( i, i ) * m_s( i ); }

            // Parte triangolare superiore
            for ( integer i = 0; i < m_nv - 1; ++i )
            {
              for ( integer j = i + 1; j < m_nv; ++j )
              {
                m_hs( i ) += m_HQ( i, j ) * m_s( j );
                m_hs( j ) += m_HQ( i, j ) * m_s( i );
              }
            }

            // Contributo da m_pq e m_xpt
            for ( integer k = 0; k < m_npt; ++k )
            {
              if ( m_pq( k ) != 0 )
              {
                Scalar temp_qts = m_pq( k ) * m_xpt.col( k ).dot( m_s );
                m_hs += temp_qts * m_xpt.col( k );
              }
            }

            // Calcola DHD, DHS, SHS (Blocco 150)
            shs = 0;
            dhs = 0;
            dhd = 0;

            for ( integer i = 0; i < m_nv; ++i )
            {
              if ( m_xbdi( i ) == 0 )
              {
                shs += m_s( i ) * m_hs( i );
                dhs += m_d( i ) * m_hs( i );
                dhd += m_d( i ) * m_hred( i );
              }
            }

            // Cerca la massima riduzione in Q (Blocco 160)
            redmax = 0;
            isav   = 0;
            redsav = 0;
            iu     = static_cast<integer>( 17.0 * angbd + 3.1 );

            for ( integer i = 1; i <= iu; ++i )
            {
              angt             = angbd * Scalar( i ) / Scalar( iu );
              sth              = ( angt + angt ) / ( 1 + angt * angt );
              Scalar temp_curv = shs + angt * ( angt * dhd - dhs - dhs );
              rednew           = sth * ( angt * dredg - sredg - Scalar( 0.5 ) * sth * temp_curv );

              if ( rednew > redmax )
              {
                redmax = rednew;
                isav   = i;
                rdprev = redsav;
              }
              else if ( i == isav + 1 ) { rdnext = rednew; }
              redsav = rednew;
            }

            // Aggiornamento D e GNEW (Blocco 170/180)
            if ( isav == 0 )
            {
              current_state = State::Final;
              break;
            }

            if ( isav < iu )
            {
              Scalar temp_interp = ( rdnext - rdprev ) / ( redmax + redmax - rdprev - rdnext );
              angt               = angbd * ( Scalar( isav ) + Scalar( 0.5 ) * temp_interp ) / Scalar( iu );
            }

            cth              = ( 1 - angt * angt ) / ( 1 + angt * angt );
            sth              = ( angt + angt ) / ( 1 + angt * angt );
            Scalar temp_curv = shs + angt * ( angt * dhd - dhs - dhs );
            sdec             = sth * ( angt * dredg - sredg - Scalar( 0.5 ) * sth * temp_curv );

            if ( sdec <= 0 )
            {
              current_state = State::Final;
              break;
            }

            // Aggiorna m_g_new e hred usando operazioni vettoriali
            m_g_new = m_g_new + ( cth - 1 ) * m_hred + sth * m_hs;
            m_hred  = cth * m_hred + sth * m_hs;

            dredg  = 0;
            gredsq = 0;

            for ( integer i = 0; i < m_nv; ++i )
            {
              if ( m_xbdi( i ) == 0 )
              {
                m_d( i ) = cth * m_d( i ) + sth * m_s( i );
                dredg += m_d( i ) * m_g_new( i );
                gredsq += m_g_new( i ) * m_g_new( i );
              }
            }

            qred += sdec;

            if ( iact > 0 && isav == iu )
            {
              // Fissa variabile, ricomincia da Boundary Loop
              ++nact;
              m_xbdi( iact - 1 ) = xsav;
              current_state      = State::Boundary_Loop;
            }
            else if ( sdec > 0.01 * qred )
            {
              // Continua il loop alternativo
              current_state = State::Alternating_Loop;
            }
            else
            {
              current_state = State::Final;
            }
            break;
          }

        // ======================================================================
        case State::Final:
          // ======================================================================
          {
            // USCITA (Blocco 190)
            m_x_new = m_x_opt + m_d;

            // Clamp ai limiti usando operazioni vettoriali
            m_x_new = m_x_new.cwiseMin( m_s_upper ).cwiseMax( m_s_lower );

            // Ricalcola d e dsq
            m_d   = m_x_new - m_x_opt;
            m_dsq = m_d.squaredNorm();

            current_state = State::Exit;
            break;
          }

        case State::Exit:
          break;
      }
    }
  }
  /**
   * @brief Updates the ZMAT and BMAT matrices after shifting an interpolation point.
   *
   * This function implements the update step of the BOBYQA algorithm when an
   * interpolation point \f$ x_k \f$ (indexed by `m_knew`) is moved to a new position.
   * It updates the factorization of the interpolation matrix and the quadratic model parameters.
   *
   * @details
   * The update procedure ensures that the new quadratic model \f$ Q_{new}(x) \f$ satisfies
   * the interpolation conditions at the new set of points. The algorithm proceeds in three phases:
   *
   * <b>Phase 1: Givens Rotations</b>
   * Transforms the matrix \f$ Z \f$ (stored in `m_zmat`) such that the row corresponding
   * to the leaving point \f$ k \f$ becomes zero everywhere except in the first column.
   * \f[
   * Z \leftarrow Z \Omega, \quad \text{s.t.} \quad Z_{k, j} = 0 \quad \forall j > 1
   * \f]
   * where \f$ \Omega \f$ is a product of Givens rotation matrices.
   *
   * <b>Phase 2: ZMAT Update</b>
   * Updates the first column of \f$ Z \f$ to account for the new point location.
   * This uses the Lagrange function values \f$ \ell_t \f$ (stored in `m_v_lag`).
   * The update formula is roughly:
   * \f[
   * z^{(1)}_{new} = \frac{\tau}{\sqrt{\sigma}} z^{(1)} - \frac{\zeta}{\sqrt{\sigma}} \ell
   * \f]
   * where \f$ \sigma \f$ is the denominator (`m_denom`) and \f$ \tau, \zeta \f$ are scalars.
   *
   * <b>Phase 3: BMAT Update</b>
   * Updates the matrix \f$ B \f$ (stored in `m_bmat`), which contains the gradient and the
   * Hessian approximation. The update minimizes the change in the Frobenius norm of the
   * Hessian approximation subject to the new interpolation conditions (Quasi-Newton update).
   * \f[
   * B_{new} = B + \Delta B
   * \f]
   * The update maintains the symmetry of the Hessian block within \f$ B \f$.
   *
   * @pre `m_knew` must be the index of the point being moved (0-based in implementation logic).
   * @pre `m_denom` (denominator \f$ \sigma \f$) must be non-zero.
   * @pre `m_zmat`, `m_bmat`, and `m_v_lag` dimensions must be consistent.
   *
   * @tparam Scalar The floating-point type (e.g., float, double).
   */
  template <typename Scalar>
  void
  BOBYQA_minimizer<Scalar>::update()
  {
    // Number of columns in ZMAT excluding the first one (which is treated specially)

    // Verify matrix dimensions in Debug mode
    assert( m_B.rows() == m_nv );
    assert( m_B.cols() == m_npt + m_nv );
    assert( m_Z.rows() == m_npt );
    assert( m_Z.cols() == m_nptm );

    // Temporary workspace vector (size: points + variables)
    Vector work( m_npt + m_nv );
    // Ensure deterministic behavior by zero-initializing the workspace
    work.setZero();

    // ========================================================================
    // PHASE 1: Givens Rotations
    // Reduces the m_knew-th row of ZMAT to a multiple of e_1.
    // ========================================================================

    // Threshold for small values to avoid unnecessary rotations (numerical stability)
    Scalar ztest = Scalar( 1e-20 ) * m_Z.cwiseAbs().maxCoeff();

    for ( integer j = 1; j < m_nptm; ++j )
    {
      Scalar const b = m_Z( m_knew, j );  // Element to annihilate

      // Apply rotation only if the element is significant
      if ( abs( b ) > ztest )
      {
        Scalar const a = m_Z( m_knew, 0 );  // Pivot element

        // Compute Givens rotation parameters: c = cos(theta), s = sin(theta)
        // r = sqrt(a^2 + b^2) computed safely via std::hypot
        Scalar const t = std::hypot( a, b );
        Scalar const c = a / t;
        Scalar const s = b / t;

        // Apply rotation to columns 0 and j:
        // [ col0' ] = [ c  s ] [ col0 ]
        // [ colj' ]   [ -s c ] [ colj ]
        Vector col0 = m_Z.col( 0 );
        Vector colj = m_Z.col( j );

        m_Z.col( 0 ) = c * col0 + s * colj;
        m_Z.col( j ) = c * colj - s * col0;

        // Explicitly zero out the target element to remove round-off error
        m_Z( m_knew, j ) = 0;
      }
    }

    // ========================================================================
    // PHASE 2: Update ZMAT (First Column)
    // Applies the rank-1 update formula to the null-space basis Z.
    // ========================================================================

    // 'a' now contains the accumulated norm of the original row due to rotations
    Scalar const a = m_Z( m_knew, 0 );

    // Compute the first part of the scalar product for the update
    work.head( m_npt ) = a * m_Z.col( 0 );
    Scalar const alpha = work( m_knew );
    Scalar const tau   = m_v_lag( m_knew );  // Lagrange multiplier for the moving point

    // Adjust VLAG temporarily for the calculation (corresponds to shifting Lagrange basis)
    m_v_lag( m_knew ) -= 1;

    Scalar const sqrt_denom = std::sqrt( m_denom );

    // Update the first column of ZMAT:
    // z_new = (tau * z - a * vlag) / sqrt(sigma)
    Scalar const tempa = tau / sqrt_denom;
    Scalar const tempb = a / sqrt_denom;
    m_Z.col( 0 )       = tempa * m_Z.col( 0 ) - tempb * m_v_lag.head( m_npt );

    // ========================================================================
    // PHASE 3: Update BMAT (Gradient and Hessian)
    // Updates BMAT using a rank-2 formula involving VLAG and WORK.
    // ========================================================================

    for ( integer j = 0; j < m_nv; ++j )
    {
      // jp is the column index corresponding to the j-th variable in the Hessian block
      integer const jp = m_npt + j;

      // Store the element B(j, knew) representing the interaction between
      // the j-th variable and the moving interpolation point
      work( jp ) = m_B( j, m_knew );

      // Calculate update coefficients based on Powell's formula
      // These mix the Lagrange vector (vlag) and the work vector
      Scalar const tempa = ( alpha * m_v_lag( jp ) - tau * work( jp ) ) / m_denom;
      Scalar const tempb = ( -m_beta * work( jp ) - tau * m_v_lag( jp ) ) / m_denom;

      m_B.row( j ).head( jp + 1 ) += tempa * m_v_lag.head( jp + 1 ) + tempb * work.head( jp + 1 );

      // Enforce symmetry in the Hessian block of BMAT.
      // If we are within the Hessian block (jp >= npt), copy the row segment
      // to the corresponding column to ensure B(j, i) == B(i, j) for the Hessian part.
      if ( jp >= m_npt )
      {
        // Copy row(j).segment(...) to col(jp).head(...)
        m_B.col( jp ).head( j + 1 ) = m_B.row( j ).segment( m_npt, j + 1 ).transpose();
      }
    }

    // Ensure the Hessian block of BMAT is symmetric (fix small inconsistencies)
    // for ( integer j = 0; j < m_nv; ++j )
    // {
    //   integer const jp = m_npt + j;
    //   for ( integer i = 0; i <= j; ++i )
    //   {
    //     Scalar a = m_B( i, jp );
    //     Scalar b = m_B( j, m_npt + i );
    //     Scalar avg = Scalar( 0.5 ) * ( a + b );
    //     m_B( i, jp ) = avg;
    //     m_B( j, m_npt + i ) = avg;
    //   }
    // }
  }

  template <typename Scalar>
  bool
  BOBYQA_minimizer<Scalar>::debug_check_bmat_symmetry( Scalar tol ) const
  {
    for ( integer j = 0; j < m_nv; ++j )
    {
      integer const jp = m_npt + j;
      for ( integer i = 0; i <= j; ++i )
      {
        if ( std::abs( m_B( i, jp ) - m_B( j, m_npt + i ) ) > tol ) return false;
      }
    }
    return true;
  }

  template <typename Scalar>
  void
  BOBYQA_minimizer<Scalar>::debug_init( integer n, integer npt )
  {
    m_nv   = n;
    m_npt  = npt;
    m_dim  = m_npt + m_nv;
    m_nptm = m_npt - m_nv - 1;

    // Resize and zero internal storage used by update
    m_B.resize( m_nv, m_npt + m_nv );
    m_B.setZero();
    m_Z.resize( m_npt, std::max<integer>( m_npt - m_nv - 1, 0 ) );
    m_Z.setZero();
    m_v_lag.resize( m_dim );
    m_v_lag.setZero();
    m_pq.resize( m_npt );
    m_pq.setZero();
    m_xpt.resize( m_nv, m_npt );
    m_xpt.setZero();
    m_HQ.resize( m_nv, m_nv );
    m_HQ.setZero();
  }

  template class BOBYQA_minimizer<double>;
  template class BOBYQA_minimizer<float>;


}  // namespace Utils
