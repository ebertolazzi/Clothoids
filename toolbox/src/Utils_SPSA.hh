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

//
// file: Utils_SPSA.hh
//

#pragma once

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#ifndef UTILS_SPSA_dot_HH
#define UTILS_SPSA_dot_HH

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wold-style-cast"
#pragma clang diagnostic ignored "-Wswitch-enum"
#pragma clang diagnostic ignored "-Wdocumentation"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
#pragma clang diagnostic ignored "-Wweak-vtables"
#pragma clang diagnostic ignored "-Wshorten-64-to-32"
#pragma clang diagnostic ignored "-Wundefined-func-template"
#pragma clang diagnostic ignored "-Wdouble-promotion"
#pragma clang diagnostic ignored "-Wsigned-enum-bitfield"
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wweak-vtables"
#pragma clang diagnostic ignored "-Wunused-template"
#pragma clang diagnostic ignored "-Wnon-virtual-dtor"
#pragma clang diagnostic ignored "-Wpadded"
#pragma clang diagnostic ignored "-Wmissing-noreturn"
#endif

#include "Utils.hh"
#include "Utils_fmt.hh"
#include "Utils_eigen.hh"

#include <algorithm>
#include <cmath>
#include <cassert>
#include <functional>
#include <limits>
#include <optional>
#include <random>
#include <utility>
#include <vector>

namespace Utils {

  using std::abs;
  using std::min;
  using std::max;

  template <typename Scalar = double>
  class SPSAGradientEstimator {
  public:
    using Vector   = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using Callback = std::function<Scalar(Vector const &, Vector *)>;

    struct Options {
      size_t   repeats            = 10;            // numero di ripetizioni (media)
      Scalar   c_base             = Scalar(1e-3);  // ampiezza base di perturbazione (se c_per_component vuoto)
      Vector   c_per_component;                    // se non vuoto, usa questi c_i (size = n)
      bool     use_rademacher     = true;          // true -> ±1, false -> gaussian N(0,1)
      bool     apply_precond      = false;         // se true, applica preconditioner diag W dopo la stima
      Vector   preconditioner;                     // vettore di pesi (size = n), usato se apply_precond == true
      unsigned int rng_seed       = std::random_device{}();
      Scalar   min_delta_abs      = Scalar(1e-12); // evita divisione per numeri troppo piccoli (gauss)
      bool     allow_gradient_arg = false;         // se true, callback può ricevere un output grad (non usato da SPSA)
    };

  private:
    Options m_opts;
    Vector  m_tmp_grad; // buffer opzionale se callback scrive gradienti (non usato per SPSA)

  public:

    SPSAGradientEstimator(Options opts = Options())
    : m_opts(std::move(opts))
    {}

    // Stima il gradiente in x. Riempie 'out_grad' (size n). Ritorna il numero di valutazioni di f effettuate.
    size_t
    estimate( Vector const & x, Callback const & f, Vector & out_grad ) {
      size_t const n = static_cast<size_t>(x.size());
      assert(n > 0);

      // prepara c_i
      Vector c;
      if ( m_opts.c_per_component.size() == 0 ) {
        c = Vector::Constant(n, m_opts.c_base);
      } else {
        assert((size_t)m_opts.c_per_component.size() == n);
        c = m_opts.c_per_component;
      }

      // preconditioner
      if (m_opts.apply_precond) {
        assert((size_t)m_opts.preconditioner.size() == n);
      }

      out_grad.setZero(n);
      if (m_opts.repeats == 0) return 0;

      // RNGs
      std::mt19937 rng(m_opts.rng_seed);
      std::uniform_int_distribution<int> rademacher_dist(0, 1); // 0 -> -1, 1 -> +1
      std::normal_distribution<Scalar> gaussian_dist(0.0, 1.0);

      Vector x_plus(n);
      Vector x_minus(n);
      Vector delta(n);
      Scalar y_plus, y_minus;

      size_t eval_count = 0;

      for ( size_t rep{0}; rep < m_opts.repeats; ++rep ) {

        // sample delta
        if (m_opts.use_rademacher) {
          for ( size_t i{0}; i < n; ++i )
            delta[i] = (rademacher_dist(rng) ? Scalar(1) : Scalar(-1));
        } else {
          for ( size_t i{0}; i < n; ++i ) {
            Scalar d = gaussian_dist(rng);
            // prevenzione di delta troppo piccolo (rare per gauss, ma utile)
            if ( abs(d) < m_opts.min_delta_abs )
              d = std::copysign( m_opts.min_delta_abs, d == 0 ? Scalar(1) : d );
            delta[i] = d;
          }
        }

        // costruisci x± = x ± c .* delta
        x_plus  = x + c.cwiseProduct(delta);
        x_minus = x - c.cwiseProduct(delta);

        // valuta funzione (notare: passiamo nullptr per gradient output; callback può ignorare secondo opzione)
        y_plus  = f(x_plus,  m_opts.allow_gradient_arg ? &m_tmp_grad : nullptr);
        y_minus = f(x_minus, m_opts.allow_gradient_arg ? &m_tmp_grad : nullptr);
        eval_count += 2;

        // aggiorna stima: g_i += (y+ - y-) / (2 * c_i * delta_i)
        const Scalar denom_factor = Scalar(0.5); // we'll multiply by (y+ - y-) * 0.5 / (c_i * delta_i)
        Scalar diff = (y_plus - y_minus);

        for (size_t i = 0; i < n; ++i) {
          Scalar di = delta[i];
          // sicurezza: se delta molto vicino a zero (solo possibile con gauss), salta o usa fallback
          if (std::abs(di) < m_opts.min_delta_abs) {
            // fallback: ignora contributo di questa ripetizione per la componente i
            continue;
          }
          out_grad[i] += diff * (denom_factor / (c[i] * di));
        }
      } // end repeats

      // media
      out_grad /= static_cast<Scalar>(m_opts.repeats);

      // applica precondizionamento (se richiesto): semplicemente scala componente-wise
      if (m_opts.apply_precond) {
        out_grad = out_grad.cwiseProduct(m_opts.preconditioner);
      }

      return eval_count;
    }

    Options const & options() const { return m_opts; }
    Options & options() { return m_opts; }
  };
}

#endif

#endif

//
// eof: Utils_LBFGS.hh
