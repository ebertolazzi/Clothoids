/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2025                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

/**
 * @file Utils_minimize_Newton.hh
 * @brief Newton optimization with Wolfe line search and adaptive regularization
 *
 * @details Implementazione di un ottimizzatore Newton regolarizzato con line search di Wolfe
 * e supporto per vincoli di scatola. Combina i vantaggi del metodo Newton
 * (convergenza quadratica locale) con la robustezza del line search per problemi non convessi.
 *
 * ## Miglioramenti Implementati
 *
 * 1. **Strategia di Regolarizzazione Robusta**:
 *    - Regola λ basandosi sulla qualità del modello quadratico
 *    - Inizializzazione adattiva basata sulla scala del problema
 *    - Gestione conservativa dei fallimenti
 *
 * 2. **Criteri di Convergenza Intelligenti**:
 *    - Separazione tra convergenza vera e stallo
 *    - Soglie relative e assolute combinate
 *    - Riconoscimento di "quasi convergenza"
 *
 * 3. **Gestione Errori Migliorata**:
 *    - Fallback graduale a gradiente quando Newton fallisce
 *    - Recovery da punti non validi
 *    - Validazione input rigorosa
 *
 * 4. **Monitoraggio Avanzato**:
 *    - Tracking del miglior punto trovato
 *    - Statistiche dettagliate delle valutazioni
 *    - Verbosità configurabile con output chiaro
 *
 * ## Strategie di Regolarizzazione
 *
 * ### Levenberg-Marquardt Adattivo:
 *    - λ regolato basandosi sul rapporto di riduzione ρ = actual_red / pred_red
 *    - ρ > 0.75 → diminuisci λ (modello accurato)
 *    - ρ < 0.1 → aumenta λ (modello inaccurato)
 *    - Altrimenti mantieni λ (modello adeguato)
 *
 * ### Inizializzazione Intelligente:
 *    - λ₀ basato sulla norma di Frobenius dell'Hessiana
 *    - Scala con la dimensione del problema
 *    - Limiti conservativi [λ_min, λ_max]
 *
 * @author Enrico Bertolazzi
 * @date 2025
 */

#pragma once

#ifndef UTILS_MINIMIZE_NEWTON_HH
#define UTILS_MINIMIZE_NEWTON_HH

#include "Utils_minimize.hh"
#include "Utils_ssolver.hh"
#include "Utils_Linesearch.hh"

namespace Utils
{

  // ===========================================================================
  // Newton Minimizer - Versione Robusta
  // ===========================================================================

  template <typename Scalar = double> class Newton_minimizer
  {
  public:
    using Vector       = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using integer      = Eigen::Index;
    using SparseMatrix = Eigen::SparseMatrix<Scalar>;
    using DenseMatrix  = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using Callback     = std::function<Scalar( Vector const &, Vector *, SparseMatrix * )>;
    using BOX          = BoxConstraintHandler<Scalar>;

    // Stati dell'ottimizzatore
    enum class Status
    {
      CONVERGED          = 0,  // Convergenza soddisfatta (‖∇f‖ < tol)
      MAX_ITERATIONS     = 1,  // Raggiunto limite iterazioni
      LINE_SEARCH_FAILED = 2,  // Line search non trova passo accettabile
      HESSIAN_FAILURE    = 3,  // Hessiana non definita positiva/singolare
      STALLED            = 4,  // Progresso insufficiente (stallo)
      INVALID_INPUT      = 5   // Input non validi (dimensioni, NaN, etc.)
    };

    // Opzioni configurabili
    struct Options
    {
      // --- Parametri di Convergenza ---
      integer max_iter = 200;    // Massimo numero di iterazioni
      Scalar  g_tol    = 1e-8;   // Tolleranza sul gradiente (norma infinito proiettata)
      Scalar  f_tol    = 1e-12;  // Tolleranza relativa sulla funzione
      Scalar  x_tol    = 1e-12;  // Tolleranza relativa sulle variabili

      // --- Regolarizzazione ---
      Scalar lambda_init = 1e-8;   // Valore iniziale di λ (0 = Newton puro)
      Scalar lambda_min  = 1e-12;  // Valore minimo di λ
      Scalar lambda_max  = 1e4;    // Valore massimo di λ
      Scalar lambda_incr = 2.0;    // Fattore di incremento (λ → λ*incr)
      Scalar lambda_decr = 0.5;    // Fattore di decremento (λ → λ*decr)

      // Parametri per strategia avanzata
      Scalar lambda_growth_factor = 5.0;   // Fattore crescita esponenziale
      Scalar lambda_decay_factor  = 0.2;   // Fattore decadimento esponenziale
      Scalar min_ratio_threshold  = 0.1;   // Soglia minima rapporto riduzione
      Scalar max_ratio_threshold  = 0.75;  // Soglia massima rapporto riduzione

      // --- Line Search ---
      StrongWolfeLineSearch<Scalar> line_search{};  // Configurazione line search

      // --- Solver ---
      bool use_dense = false;  // Usa solver denso invece che sparso

      // --- Verbosità ---
      integer verbosity = 1;  // 0: silenzioso, 1: base, 2: dettagliato, 3: debug

      // --- Strategie Avanzate ---
      bool   adaptive_lambda_init     = true;  // Inizializza λ basandosi su ‖H‖
      Scalar lambda_init_factor       = 1e-6;  // Fattore per λ_init = factor * ‖H‖_F
      bool   use_gradient_fallback    = true;  // Fallback a gradiente se Newton fallisce
      bool   keep_best_point          = true;  // Mantieni il miglior punto trovato
      Scalar quasi_convergence_factor = 10.0;  // Fattore per "quasi convergenza"
    };

  private:
    Options m_opts;                                          // Opzioni configurabili
    BOX     m_box;                                           // Gestore vincoli di scatola
    Scalar  m_eps = std::numeric_limits<Scalar>::epsilon();  // Precisione di macchina

    // --- Stato Corrente ---
    Scalar       m_f;       // Valore funzione corrente
    Scalar       m_lambda;  // Parametro di regolarizzazione corrente
    Scalar       m_gnorm;   // Norma del gradiente proiettato
    Vector       m_x;       // Punto corrente
    Vector       m_g;       // Gradiente corrente
    SparseMatrix m_H;       // Hessiana corrente (sparsa)

    // --- Miglior Punto Trovato ---
    Vector m_best_x;                                       // Miglior punto trovato
    Scalar m_best_f = std::numeric_limits<Scalar>::max();  // Miglior valore funzione trovato

    // --- Statistiche ---
    Status  m_status       = Status::STALLED;  // Stato finale
    integer m_iter         = 0;                // Iterazioni effettuate
    integer m_f_evals      = 0;                // Valutazioni funzione
    integer m_h_evals      = 0;                // Valutazioni Hessiana
    integer m_ls_evals     = 0;                // Valutazioni line search
    integer m_newton_fails = 0;                // Fallimenti direzione Newton
    Scalar  m_f_init;                          // Valore iniziale funzione

    // --- Metriche di Progresso ---
    Scalar m_last_df;                                          // Ultima riduzione funzione
    Scalar m_last_dx;                                          // Ultima variazione variabili
    Scalar m_best_gnorm = std::numeric_limits<Scalar>::max();  // Migliore norma gradiente raggiunta

    // =========================================================================
    // METODI PRIVATI - IMPLEMENTAZIONE DETTAGLIATA
    // =========================================================================

    /**
     * @brief Inizializza λ basandosi sulla norma dell'Hessiana
     *
     * Strategia: λ₀ = ε * ‖H‖_F / (1 + ‖x‖) dove ε ≈ 1e-6
     * Giustificazione: regolarizzazione proporzionale alla scala dell'Hessiana
     * e inversamente proporzionale alla scala delle variabili
     *
     * @param H Hessiana iniziale
     * @return Scalar Valore iniziale per λ
     */
    Scalar compute_initial_lambda( SparseMatrix const & H ) const
    {
      if ( !m_opts.adaptive_lambda_init ) return m_opts.lambda_init;

      // Calcola norma di Frobenius dell'Hessiana
      Scalar norm_H = 0;
      for ( integer j = 0; j < H.cols(); ++j )
      {
        for ( typename SparseMatrix::InnerIterator it( H, j ); it; ++it ) { norm_H += it.value() * it.value(); }
      }
      norm_H = std::sqrt( norm_H );

      // Scala con la norma delle variabili
      Scalar norm_x = m_x.norm();

      // λ₀ = fattore * norma_Hessiana / (1 + norma_variabili)
      Scalar lambda = m_opts.lambda_init_factor * norm_H / ( 1.0 + norm_x );

      // Limita tra min e max
      lambda = std::max( m_opts.lambda_min, std::min( lambda, m_opts.lambda_max ) );

      if ( m_opts.verbosity >= 2 )
        fmt::print( PrintColors::INFO, "    λ iniziale: {:.2e} (‖H‖={:.2e}, ‖x‖={:.2e})\n", lambda, norm_H, norm_x );

      return lambda;
    }

    /**
     * @brief Calcola direzione di Newton regolarizzata (H + λI)p = -g
     *
     * Implementa due strategie:
     * 1. Solver denso (use_dense = true): converte in densa, usa decomposizione LDLT
     * 2. Solver sparso (default): usa solver simmetrico per matrici sparse
     *
     * Controlli:
     * - Finitezza della soluzione (no NaN/Inf)
     * - Direzione di discesa (p·g < 0)
     *
     * @param dir [out] Direzione calcolata
     * @return true se successo, false altrimenti
     */
    bool compute_direction( Vector & dir )
    {
      try
      {
        // CASO 1: SOLVER DENSO (per piccoli problemi o matrici dense)
        if ( m_opts.use_dense )
        {
          DenseMatrix                  H_dense = m_H;  // Conversione a densa
          DenseSymmetricSolver<Scalar> solver( H_dense, m_lambda );

          // Risolvi (H + λI)p = -g
          dir.noalias() = -solver.solve( m_g );

          // Controlla validità numerica
          if ( !dir.allFinite() )
          {
            if ( m_opts.verbosity >= 3 ) fmt::print( PrintColors::WARNING, "    Solver denso prodotto NaN/Inf\n" );
            return false;
          }
        }
        // CASO 2: SOLVER SPARSO (default per grandi problemi)
        else
        {
          SparseSymmetricSolver<Scalar> solver( m_H, m_lambda );
          dir.noalias() = -solver.solve( m_g );

          if ( !dir.allFinite() )
          {
            if ( m_opts.verbosity >= 3 ) fmt::print( PrintColors::WARNING, "    Solver sparso prodotto NaN/Inf\n" );
            return false;
          }
        }

        // VERIFICA DIREZIONE DI DISCESA
        // La direzione deve soddisfare p·g < 0 (discesa)
        Scalar slope     = dir.dot( m_g );
        Scalar slope_tol = -m_eps * ( 1.0 + std::abs( slope ) );

        if ( slope > slope_tol )
        {
          // Direzione non di discesa
          if ( m_opts.verbosity >= 3 )
          {
            fmt::print(
              PrintColors::WARNING,
              "    Direzione non di discesa: p·g = {:.2e} > {:.2e}\n",
              slope,
              slope_tol );
          }
          return false;
        }

        return true;
      }
      catch ( std::exception const & e )
      {
        // Eccezione dal solver (matrice singolare, etc.)
        if ( m_opts.verbosity >= 3 ) fmt::print( PrintColors::WARNING, "    Eccezione solver: {}\n", e.what() );
        return false;
      }
      catch ( ... )
      {
        if ( m_opts.verbosity >= 3 ) fmt::print( PrintColors::WARNING, "    Eccezione sconosciuta nel solver\n" );
        return false;
      }
    }

    /**
     * @brief Strategia di regolarizzazione con tentativi multipli e fallback graduale
     *
     * Tentativi in sequenza con λ crescenti:
     * 1. λ corrente
     * 2. λ * incr
     * 3. λ * incr²
     * 4. λ_max
     *
     * Se tutti falliscono, fallback a gradiente (se configurato)
     *
     * @param dir [out] Direzione calcolata
     * @return true se trovata direzione valida, false altrimenti
     */
    bool try_regularized_direction( Vector & dir )
    {
      // Primo tentativo con λ corrente
      if ( compute_direction( dir ) ) return true;

      // Sequenza di λ crescenti (crescita esponenziale moderata)
      std::array<Scalar, 4> lambdas = {
        m_lambda * m_opts.lambda_incr,
        m_lambda * m_opts.lambda_incr * m_opts.lambda_incr,
        m_opts.lambda_max,
        m_opts.lambda_max  // Ultimo tentativo con λ massimo
      };

      for ( size_t i = 0; i < lambdas.size(); ++i )
      {
        m_lambda = lambdas[i];

        if ( m_opts.verbosity >= 3 ) fmt::print( "      Tentativo direzione con λ = {:.2e}\n", m_lambda );

        if ( compute_direction( dir ) )
        {
          if ( m_opts.verbosity >= 2 )
          {
            fmt::print( PrintColors::INFO, "    Direzione accettata con λ = {:.2e} (tentativo {})\n", m_lambda, i + 2 );
          }
          return true;
        }
      }

      // Tutti i tentativi falliti, fallback a gradiente se configurato
      if ( m_opts.use_gradient_fallback )
      {
        dir = -m_g;
        m_box.project_direction( m_x, dir );
        ++m_newton_fails;

        // Verifica che sia una direzione di discesa
        Scalar slope = dir.dot( m_g );
        if ( slope < -m_eps * ( 1.0 + std::abs( slope ) ) )
        {
          if ( m_opts.verbosity >= 2 ) fmt::print( PrintColors::WARNING, "    Fallback a direzione di gradiente\n" );
          return true;
        }
      }

      if ( m_opts.verbosity >= 1 )
        fmt::print( PrintColors::WARNING, "  Fallimento direzione con tutti i λ e gradiente\n" );

      return false;
    }

    /**
     * @brief Calcola norma del gradiente proiettato (norma infinito)
     *
     * Per problemi vincolati, usa la norma infinito del gradiente proiettato:
     * ‖∇f‖_∞,proj = max_i |∂f/∂x_i| con x_i al bordo se attivo
     *
     * Per problemi non vincolati: ‖∇f‖_2
     *
     * @return Scalar Norma del gradiente proiettato
     */
    Scalar projected_grad_norm() const
    {
      if ( m_box.is_active() )
        return m_box.projected_gradient_norm_inf( m_x, m_g );
      else
        return m_g.template lpNorm<Eigen::Infinity>();
    }

    /**
     * @brief Aggiorna λ basandosi sulla qualità del modello quadratico
     *
     * Strategia basata su Levenberg-Marquardt/Trust Region:
     * ρ = riduzione_effettiva / riduzione_predetta
     *
     * Classificazione:
     * - ρ > 0.75: modello eccellente → diminuisci λ
     * - ρ ∈ [0.1, 0.75]: modello adeguato → mantieni λ
     * - ρ < 0.1: modello scarso → aumenta λ
     * - ρ ≤ 0: modello catastrofico → aumenta λ
     *
     * @param actual_red Riduzione effettiva f(x) - f(x+αp)
     * @param alpha Lunghezza passo
     * @param p Direzione di ricerca
     */
    void adjust_lambda( Scalar actual_red, Scalar alpha, Vector const & p )
    {
      // Calcola riduzione predetta dal modello quadratico:
      // pred_red = -αgᵀp - ½α²pᵀHp
      Scalar gTp      = m_g.dot( p );
      Scalar pHp      = p.dot( m_H * p );
      Scalar pred_red = -alpha * ( gTp + 0.5 * alpha * pHp );

      // Caso patologico: riduzione predetta non positiva
      if ( pred_red <= 0 )
      {
        // Modello cattivo → aumenta regolarizzazione moderatamente
        m_lambda = std::min( m_opts.lambda_max, m_opts.lambda_growth_factor * m_lambda );
        if ( m_opts.verbosity >= 3 ) fmt::print( "      λ↑ (pred_red ≤ 0): {:.2e}\n", m_lambda );
        return;
      }

      // Rapporto riduzione effettiva/predetta
      Scalar ratio = actual_red / pred_red;

      // Regola λ in modo conservativo
      if ( ratio > m_opts.max_ratio_threshold )
      {
        // Modello eccellente → diminuisci λ gradualmente
        m_lambda = std::max( m_opts.lambda_min, m_opts.lambda_decay_factor * m_lambda );
        if ( m_opts.verbosity >= 3 ) fmt::print( "      λ↓ (ratio={:.2f}): {:.2e}\n", ratio, m_lambda );
      }
      else if ( ratio < m_opts.min_ratio_threshold )
      {
        // Modello scarso → aumenta λ moderatamente
        m_lambda = std::min( m_opts.lambda_max, m_opts.lambda_growth_factor * m_lambda );
        if ( m_opts.verbosity >= 3 ) fmt::print( "      λ↑ (ratio={:.2f}): {:.2e}\n", ratio, m_lambda );
      }
      else
      {
        // Modello adeguato → mantieni λ
        if ( m_opts.verbosity >= 3 ) fmt::print( "      λ mantiene (ratio={:.2f})\n", ratio );
      }
    }

    /**
     * @brief Controlla se l'ottimizzazione è convergenza
     *
     * Criteri (in ordine di priorità):
     * 1. ‖∇f‖_proj ≤ g_tol (CONVERGENZA VERA)
     * 2. ‖∇f‖_proj ≤ quasi_convergence_factor * g_tol E (Δf_rel ≤ f_tol O Δx_rel ≤ x_tol)
     * 3. Δf_rel ≤ f_tol E Δx_rel ≤ x_tol E ‖∇f‖_proj < 1e-4 (STALLO)
     *
     * @param df Variazione assoluta funzione
     * @param dx Variazione assoluta variabili (norma passo)
     * @return true se convergenza/stallo, false altrimenti
     */
    bool converged( Scalar df, Scalar dx ) const
    {
      // 1. GRADIENTE SUFFICIENTEMENTE PICCOLO (CONVERGENZA VERA)
      if ( m_gnorm <= m_opts.g_tol )
      {
        if ( m_opts.verbosity >= 2 )
          fmt::print( PrintColors::SUCCESS, "    ‖∇f‖ = {:.2e} ≤ {:.2e}\n", m_gnorm, m_opts.g_tol );
        return true;
      }

      // 2. QUASI CONVERGENZA (gradiente moderatamente piccolo e progresso insufficiente)
      if ( m_gnorm <= m_opts.quasi_convergence_factor * m_opts.g_tol )
      {
        Scalar f_rel = std::abs( df ) / ( 1.0 + std::abs( m_f ) );
        Scalar x_rel = dx / ( 1.0 + m_x.norm() );

        if ( f_rel < m_opts.f_tol || x_rel < m_opts.x_tol )
        {
          if ( m_opts.verbosity >= 2 )
          {
            fmt::print(
              PrintColors::INFO,
              "    Quasi convergenza: ‖∇f‖={:.2e}, Δf_rel={:.2e}, Δx_rel={:.2e}\n",
              m_gnorm,
              f_rel,
              x_rel );
          }
          return true;
        }
      }

      // 3. STALLO VERO (progresso infinitesimo e gradiente non troppo grande)
      Scalar f_rel = std::abs( df ) / ( 1.0 + std::abs( m_f ) );
      Scalar x_rel = dx / ( 1.0 + m_x.norm() );

      if ( f_rel < m_opts.f_tol && x_rel < m_opts.x_tol && m_gnorm < 1e-4 )
      {
        if ( m_opts.verbosity >= 2 )
        {
          fmt::print(
            PrintColors::WARNING,
            "    Stallo: ‖∇f‖={:.2e}, Δf_rel={:.2e}, Δx_rel={:.2e}\n",
            m_gnorm,
            f_rel,
            x_rel );
        }
        return true;
      }

      return false;
    }

    /**
     * @brief Determina lo stato finale dell'ottimizzazione
     *
     * Analizza le condizioni di terminazione per assegnare lo status:
     * - CONVERGED: gradiente sotto tolleranza
     * - STALLED: gradiente sopra tolleranza ma progresso insufficiente
     * - MAX_ITERATIONS: superato limite iterazioni
     * - Altri: fallimenti specifici
     *
     * @return Status appropriato
     */
    Status determine_final_status() const
    {
      // Convergenza vera (gradiente piccolo)
      if ( m_gnorm <= m_opts.g_tol ) return Status::CONVERGED;

      // Quasi convergenza (gradiente moderatamente piccolo e progresso insufficiente)
      // Questa è una forma di convergenza accettabile
      if ( m_gnorm <= m_opts.quasi_convergence_factor * m_opts.g_tol ) return Status::CONVERGED;

      // Limite iterazioni raggiunto
      if ( m_iter >= m_opts.max_iter ) return Status::MAX_ITERATIONS;

      // Stallo (progresso insufficiente ma non convergenza)
      return Status::STALLED;
    }

    // =========================================================================
    // METODI DI UTILITY PER OUTPUT
    // =========================================================================

    void print_iter( integer it, Scalar f, Scalar gn, Scalar a, Scalar df ) const
    {
      if ( m_opts.verbosity < 2 ) return;

      // Colori basati sul progresso
      auto color = PrintColors::INFO;
      if ( df > 0.1 * std::abs( m_f_init ) )
        color = PrintColors::SUCCESS;
      else if ( df < 0 )
        color = PrintColors::WARNING;

      fmt::print(
        color,
        "[{:4d}] F={:12.6e} ΔF={:9.2e} ‖g‖={:8.2e} λ={:8.2e} α={:8.2e}\n",
        it,
        f,
        df,
        gn,
        m_lambda,
        a );
    }

    void print_header( integer n ) const
    {
      if ( m_opts.verbosity < 1 ) return;

      fmt::print( "╔══════════════════════════════════════════════════════════════════════╗\n" );
      fmt::print( "║              NEWTON OPTIMIZER WITH ANALYTICAL HESSIAN               ║\n" );
      fmt::print( "╠══════════════════════════════════════════════════════════════════════╣\n" );
      fmt::print(
        "║ Dimension: {:4d}        MaxIter: {:4d}        GTol: {:8.2e}       ║\n",
        n,
        m_opts.max_iter,
        m_opts.g_tol );
      fmt::print(
        "║ λ₀={:<8.2e}     Range=[{:<8.2e}, {:<8.2e}]                      ║\n",
        m_opts.lambda_init,
        m_opts.lambda_min,
        m_opts.lambda_max );
      if ( m_box.is_active() )
      {
        integer n_active = m_box.num_active( m_x );
        fmt::print( "║ Box constraints: ACTIVE ({:d}/{:d} variabili vincolate)              ║\n", n_active, n );
      }
      else
        fmt::print( "║ Box constraints: INACTIVE                                            ║\n" );
      fmt::print( "║ Initial F = {:<12.6e}                                             ║\n", m_f_init );
      fmt::print( "╚══════════════════════════════════════════════════════════════════════╝\n" );
      if ( m_opts.verbosity >= 2 ) fmt::print( " Iter     F(x)          ΔF         ‖∇f‖         λ          α\n" );
    }

    void print_summary() const
    {
      if ( m_opts.verbosity < 1 ) return;

      // Colore basato sullo stato
      auto color = PrintColors::INFO;
      if ( m_status == Status::CONVERGED )
        color = PrintColors::SUCCESS;
      else if ( m_status == Status::STALLED || m_status == Status::MAX_ITERATIONS )
        color = PrintColors::WARNING;
      else
        color = PrintColors::ERROR;

      fmt::print( "\n╔══════════════════════════════════════════════════════════════════════╗\n" );
      fmt::print( "║                     OPTIMIZATION SUMMARY                            ║\n" );
      fmt::print( "╠══════════════════════════════════════════════════════════════════════╣\n" );
      fmt::print( color, "║ Status:            {:<45} ║\n", to_string( m_status ) );
      fmt::print( "║ Iterations:        {:<45d} ║\n", m_iter );
      fmt::print( "║ Function evals:    {:<45d} ║\n", m_f_evals );
      fmt::print( "║ Hessian evals:     {:<45d} ║\n", m_h_evals );
      fmt::print( "║ Line search evals: {:<45d} ║\n", m_ls_evals );
      fmt::print( "║ Newton fails:      {:<45d} ║\n", m_newton_fails );
      fmt::print( "╠══════════════════════════════════════════════════════════════════════╣\n" );
      fmt::print( "║ Final ‖∇f‖:        {:<45.2e} ║\n", m_gnorm );
      fmt::print( "║ Final f:           {:<45.6e} ║\n", m_f );
      fmt::print( "║ Initial f:         {:<45.6e} ║\n", m_f_init );
      fmt::print( "║ Reduction:         {:<45.6e} ║\n", m_f_init - m_f );
      if ( m_opts.keep_best_point && m_best_f < m_f - m_eps )
        fmt::print( "║ Best f found:      {:<45.6e} ║\n", m_best_f );
      fmt::print( "╚══════════════════════════════════════════════════════════════════════╝\n" );
    }

  public:
    // =========================================================================
    // COSTRUTTORE E METODI PUBBLICI
    // =========================================================================

    /**
     * @brief Costruttore
     * @param opts Opzioni di configurazione
     */
    Newton_minimizer( Options opts = Options() ) : m_opts( opts ), m_lambda( opts.lambda_init ) {}

    /**
     * @brief Imposta vincoli di scatola
     * @param lower Limiti inferiori
     * @param upper Limiti superiori
     */
    void set_bounds( Vector const & lower, Vector const & upper ) { m_box.set_bounds( lower, upper ); }

    /**
     * @brief Converti Status in stringa
     */
    static std::string to_string( Status s )
    {
      switch ( s )
      {
        case Status::CONVERGED: return "CONVERGED";
        case Status::MAX_ITERATIONS: return "MAX_ITERATIONS";
        case Status::LINE_SEARCH_FAILED: return "LINE_SEARCH_FAILED";
        case Status::HESSIAN_FAILURE: return "HESSIAN_FAILURE";
        case Status::STALLED: return "STALLED";
        case Status::INVALID_INPUT: return "INVALID_INPUT";
        default: return "UNKNOWN";
      }
    }

    /**
     * @brief Esegue l'ottimizzazione
     *
     * Algoritmo principale:
     * 1. Inizializzazione e valutazione iniziale
     * 2. Loop principale (fino a max_iter o convergenza)
     *    a. Calcolo direzione di Newton regolarizzata
     *    b. Line search di Wolfe
     *    c. Aggiornamento punto e parametri
     * 3. Finalizzazione e output
     *
     * @param x0 Punto iniziale
     * @param callback Funzione che valuta f, ∇f, ∇²f
     */
    void minimize( Vector const & x0, Callback const & callback )
    {
      // === INIZIALIZZAZIONE ===
      m_iter         = 0;
      m_f_evals      = 0;
      m_h_evals      = 0;
      m_ls_evals     = 0;
      m_newton_fails = 0;
      m_x            = x0;
      m_g.resize( x0.size() );
      m_H.resize( x0.size(), x0.size() );

      // Validazione input
      if ( !m_x.allFinite() )
      {
        m_status = Status::INVALID_INPUT;
        if ( m_opts.verbosity >= 1 ) fmt::print( PrintColors::ERROR, "ERRORE: Punto iniziale contiene NaN/Inf\n" );
        return;
      }

      // Proiezione sui vincoli
      m_box.project( m_x );

      // Valutazione iniziale
      m_f       = callback( m_x, &m_g, &m_H );
      m_f_evals = m_h_evals = 1;

      if ( !std::isfinite( m_f ) || !m_g.allFinite() )
      {
        m_status = Status::INVALID_INPUT;
        if ( m_opts.verbosity >= 1 ) fmt::print( PrintColors::ERROR, "ERRORE: Valori non finiti in f/g iniziali\n" );
        return;
      }

      // Inizializzazione stato
      m_best_f = m_f_init = m_f;
      m_best_x            = m_x;
      m_gnorm             = projected_grad_norm();
      m_best_gnorm        = m_gnorm;

      // Inizializzazione λ adattativa
      if ( m_opts.adaptive_lambda_init ) m_lambda = compute_initial_lambda( m_H );

      // Reset dello stato a STALLED (in esecuzione)
      m_status = Status::STALLED;

      print_header( x0.size() );

      // === LOOP PRINCIPALE ===
      for ( m_iter = 1; m_iter <= m_opts.max_iter; ++m_iter )
      {
        // Controllo convergenza primario (gradiente)
        if ( m_gnorm <= m_opts.g_tol )
        {
          m_status = Status::CONVERGED;
          if ( m_opts.verbosity >= 1 )
            fmt::print(
              PrintColors::SUCCESS,
              "\nConvergenza raggiunta: ‖∇f‖ = {:.2e} ≤ {:.2e}\n",
              m_gnorm,
              m_opts.g_tol );
          break;
        }

        // --- CALCOLO DIREZIONE ---
        Vector p;
        if ( !try_regularized_direction( p ) )
        {
          m_status = Status::HESSIAN_FAILURE;
          if ( m_opts.verbosity >= 1 ) fmt::print( PrintColors::ERROR, "\nFallimento direzione Newton\n" );
          break;
        }

        // Proiezione direzione sui vincoli attivi
        m_box.project_direction( m_x, p );

        // Verifica finale direzione di discesa
        Scalar g0p       = p.dot( m_g );
        Scalar min_slope = -m_eps * ( 1.0 + std::abs( g0p ) );

        if ( g0p >= min_slope )
        {
          // Direzione non di discesa sufficiente → fallback a gradiente
          p = -m_g;
          m_box.project_direction( m_x, p );
          g0p = p.dot( m_g );
          ++m_newton_fails;

          if ( m_opts.verbosity >= 2 )
            fmt::print( PrintColors::WARNING, "    Fallback a gradiente (g·p={:.2e})\n", g0p );
        }

        // --- LINE SEARCH ---
        auto ls_result = m_opts.line_search(
          m_f,
          g0p,
          m_x,
          p,
          [&]( Vector const & x_eval, Vector * g_eval = nullptr ) -> Scalar
          {
            // Funzione wrapper per valutazioni line search
            Vector x_proj = x_eval;
            m_box.project( x_proj );

            // Valuta funzione (e opzionalmente gradiente)
            Scalar f_val;
            if ( g_eval != nullptr )
              f_val = callback( x_proj, g_eval, nullptr );
            else
              f_val = callback( x_proj, nullptr, nullptr );

            ++m_ls_evals;
            return f_val;
          } );

        if ( !ls_result )
        {
          m_status = Status::LINE_SEARCH_FAILED;
          if ( m_opts.verbosity >= 1 ) fmt::print( PrintColors::ERROR, "\nLine search fallito\n" );
          break;
        }

        // Estrai risultati line search
        auto [alpha, ls_cnt] = *ls_result;
        m_f_evals += ls_cnt;

        // --- AGGIORNAMENTO PUNTO ---
        Vector x_new = m_x + alpha * p;
        m_box.project( x_new );

        Vector g_new( m_x.size() );
        Scalar f_new = callback( x_new, &g_new, &m_H );
        ++m_f_evals;
        ++m_h_evals;

        // Verifica validità nuovo punto
        if ( !std::isfinite( f_new ) || !g_new.allFinite() )
        {
          // Punto non valido → aumenta regolarizzazione e continua
          m_lambda = std::min( m_opts.lambda_max, m_opts.lambda_growth_factor * m_lambda );
          if ( m_opts.verbosity >= 2 )
            fmt::print( PrintColors::WARNING, "    Punto non valido, λ↑ = {:.2e}\n", m_lambda );
          continue;
        }

        // Calcola progresso
        Scalar df = m_f - f_new;
        Scalar dx = ( alpha * p ).norm();

        // Accetta step se migliora (con tolleranza numerica)
        if ( f_new < m_f - m_eps * ( 1.0 + std::abs( m_f ) ) )
        {
          // Aggiorna stato
          m_f       = f_new;
          m_x       = x_new;
          m_g       = g_new;
          m_last_df = df;
          m_last_dx = dx;

          // Aggiorna miglior punto trovato
          if ( f_new < m_best_f )
          {
            m_best_f = f_new;
            m_best_x = x_new;
          }

          // Aggiorna norma gradiente
          m_gnorm = projected_grad_norm();
          if ( m_gnorm < m_best_gnorm ) m_best_gnorm = m_gnorm;

          // Regola λ basandosi sulla qualità dello step
          adjust_lambda( df, alpha, p );

          // Output iterazione
          print_iter( m_iter, m_f, m_gnorm, alpha, df );

          // Controlla convergenza/stallo
          if ( converged( df, dx ) )
          {
            m_status = determine_final_status();
            break;
          }
        }
        else
        {
          // Step rifiutato (non sufficiente riduzione)
          // Aumenta regolarizzazione per step più conservativi
          m_lambda = std::min( m_opts.lambda_max, m_opts.lambda_growth_factor * m_lambda );

          if ( m_opts.verbosity >= 2 )
            fmt::print( PrintColors::WARNING, "    Step rifiutato (Δf={:.2e}), λ↑ = {:.2e}\n", df, m_lambda );
        }
      }

      // === FINALIZZAZIONE ===

      // Controllo limite iterazioni - SOLO se non abbiamo già uno stato determinato
      if ( m_status == Status::STALLED && m_iter > m_opts.max_iter ) { m_status = Status::MAX_ITERATIONS; }

      // Ripristina miglior punto trovato (se configurato)
      if ( m_opts.keep_best_point && m_best_f < m_f - m_eps )
      {
        m_x = m_best_x;
        m_f = m_best_f;
        // Rivaluta gradiente per consistenza
        m_f = callback( m_x, &m_g, &m_H );
        ++m_f_evals;
        ++m_h_evals;
        m_gnorm = projected_grad_norm();
      }

      // Output riepilogativo
      print_summary();
    }

    // =========================================================================
    // METODI DI ACCESSO
    // =========================================================================

    Status         status() const { return m_status; }
    integer        iterations() const { return m_iter; }
    integer        function_evals() const { return m_f_evals; }
    integer        hessian_evals() const { return m_h_evals; }
    integer        line_search_evals() const { return m_ls_evals; }
    integer        newton_failures() const { return m_newton_fails; }
    Scalar         final_f() const { return m_f; }
    Scalar         initial_f() const { return m_f_init; }
    Scalar         final_grad_norm() const { return m_gnorm; }
    Scalar         best_grad_norm() const { return m_best_gnorm; }
    Vector const & solution() const { return m_best_x; }
    Scalar         current_lambda() const { return m_lambda; }

    void use_dense_solver( bool use_dense ) { m_opts.use_dense = use_dense; }
    void set_verbosity( integer level ) { m_opts.verbosity = level; }
    void set_lambda( Scalar lambda ) { m_lambda = lambda; }

    /**
     * @brief Resetta lo stato dell'ottimizzatore
     */
    void reset()
    {
      m_status = Status::STALLED;
      m_iter = m_f_evals = m_h_evals = m_ls_evals = m_newton_fails = 0;
      m_lambda                                                     = m_opts.lambda_init;
      m_best_f                                                     = std::numeric_limits<Scalar>::max();
      m_best_gnorm                                                 = std::numeric_limits<Scalar>::max();
    }
  };

}  // namespace Utils

#endif  // UTILS_MINIMIZE_NEWTON_HH
