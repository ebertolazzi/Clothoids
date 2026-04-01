/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                      |
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
// file: Utils_search_intervals.hh
//

#pragma once

#ifndef UTILS_SEARCH_INTERVALS_HH
#define UTILS_SEARCH_INTERVALS_HH

#include "Utils.hh"
#include <algorithm>  // per std::upper_bound

namespace Utils
{

  // ============================================================================
  // INTERVAL SEARCH FUNCTIONS
  // ============================================================================

  /**
   * \brief Searches for the interval containing a value in a sorted array
   * \tparam T_int Integer type for indices
   * \tparam T_real Real type for array values
   * \param npts Number of points in the array (must be >= 2)
   * \param X Sorted array of real numbers
   * \param x Value to locate (may be adjusted if closed range)
   * \param last_interval Reference to last known interval index, updated with found interval
   * \param closed If true, treats array as closed periodic range
   * \param can_extend If false, x must be within array bounds
   */
  template <typename T_int, typename T_real> inline void search_interval(
    T_int        npts,
    T_real const X[],
    T_real &     x,
    T_int &      last_interval,
    bool         closed,
    bool         can_extend )
  {
    using std::upper_bound;

    // Validazione iniziale
    T_int n{ npts - 1 };
    UTILS_ASSERT(
      npts > 1 && last_interval >= 0 && last_interval < n,
      "In search_interval( npts={}, X, x={}, last_interval={}, "
      "closed={}, can_extend={})\n"
      "npts must be >= 2 and last_interval must be in [0,npts-2]\n",
      npts,
      x,
      last_interval,
      closed,
      can_extend );

    T_real xl{ X[0] };
    T_real xr{ X[n] };

    if ( closed )
    {
      T_real L{ xr - xl };

      // Caso degenere: intervallo di lunghezza zero
      if ( L <= 0 )
      {
        last_interval = 0;
        return;
      }

      // Normalizza x nell'intervallo [xl, xr)
      T_real t{ x - xl };
      t = fmod( t, L );
      if ( t < 0 ) t += L;

      // CORREZIONE CRITICA: Gestione precisa dei bordi
      // Dopo il wrapping, t è in [0, L)
      // Se t == 0, x è esattamente sul bordo sinistro (xl) o un suo multiplo
      // Per coerenza ciclica, deve appartenere al PRIMO intervallo [0,1]
      // a meno che x non sia esattamente xr (nel qual caso è equivalente a xl)

      x = xl + t;

      // Gestione speciale per x esattamente uguale a xr
      // (dopo il wrapping, t=0 ma il valore originale era xr)
      // Lo trattiamo come se fosse xl, che va nel primo intervallo
      // MA: c'è un caso speciale: se x è esattamente xr, dobbiamo decidere
      // se assegnarlo all'ultimo intervallo o al primo.
      // Per coerenza con la classe, probabilmente va al primo intervallo.
    }
    else
    {
      UTILS_ASSERT(
        can_extend || ( x >= xl && x <= xr ),
        "In search_interval( npts={}, X, x={}, last_interval={}, "
        "closed={}, can_extend={})\n"
        "out of range: [{},{}]\n",
        npts,
        x,
        last_interval,
        closed,
        can_extend,
        xl,
        xr );
    }

    // === OTTIMIZZAZIONE: Ricerca adattiva basata sulla posizione corrente ===
    T_real const * XL{ X + last_interval };

    // Caso 1: x è a destra dell'intervallo corrente
    if ( x >= XL[1] )
    {
      if ( last_interval < n - 1 )
      {
        // Check rapido: intervallo successivo?
        if ( x < XL[2] ) { ++last_interval; }
        else
        {
          // Binary search nella parte destra
          T_real const * XE = X + n + 1;
          last_interval     = T_int( upper_bound( XL + 2, XE, x ) - X ) - 1;

          // Clamp per sicurezza
          if ( last_interval > n - 1 ) last_interval = n - 1;
        }
      }
      // else: già all'ultimo intervallo, rimani lì
    }
    // Caso 2: x è a sinistra dell'intervallo corrente
    else if ( x < XL[0] )
    {
      // Check rapido: intervallo precedente?
      if ( last_interval > 0 && x >= XL[-1] ) { --last_interval; }
      else
      {
        // Binary search nella parte sinistra
        last_interval = T_int( upper_bound( X, XL, x ) - X ) - 1;

        // Clamp per sicurezza
        if ( last_interval < 0 ) last_interval = 0;
      }
    }
    // Caso 3: x è già nell'intervallo corretto [XL[0], XL[1])
    // => Non fare nulla

    // === OTTIMIZZAZIONE: Protezione intervalli degeneri ===
    // Cerca il primo intervallo non degenere a sinistra se necessario
    if ( last_interval < n && X[last_interval + 1] <= X[last_interval] )
    {
      // L'intervallo è degenere o ha larghezza zero
      // Cerca il primo intervallo valido a sinistra
      T_int i = last_interval;
      while ( i > 0 && X[i] <= X[i - 1] ) --i;
      last_interval = i;
    }

    // CORREZIONE FINALE: Gestione del bordo destro in modalità chiusa
    if ( closed )
    {
      // Se x è esattamente uguale a xr (dopo il wrapping è xl)
      // e siamo nell'ultimo intervallo, dobbiamo spostarci al primo intervallo
      // perché xl appartiene al primo intervallo [0,1] per coerenza ciclica
      if ( x == xl && last_interval == n - 1 ) { last_interval = 0; }
    }

    // Validazione finale
    UTILS_ASSERT(
      last_interval >= 0 && last_interval < n,
      "In search_interval( npts={}, X, x={}, last_interval={}, "
      "closed={}, can_extend={})\n"
      "computed last_interval={} out of range [0,{}]\n",
      npts,
      x,
      last_interval,
      closed,
      can_extend,
      last_interval,
      n - 1 );
  }
  /**
   * \brief Legacy wrapper (deprecated)
   */
  template <typename T_int, typename T_real> inline void searchInterval(
    T_int        npts,
    T_real const X[],
    T_real &     x,
    T_int &      last_interval,
    bool         closed,
    bool         can_extend )
  {
    search_interval( npts, X, x, last_interval, closed, can_extend );
  }

}  // namespace Utils

#endif  // UTILS_SEARCH_INTERVALS_HH
