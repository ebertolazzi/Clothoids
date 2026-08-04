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
// file: Utils_progress_bar.hh
//

#pragma once

#ifndef UTILS_PROGRESS_BAR_HH
#define UTILS_PROGRESS_BAR_HH

#include <algorithm>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include "Utils.hh"
#include "Utils_fmt.hh"

namespace Utils
{

  using std::ceil;
  using std::floor;
  using std::round;

  namespace detail
  {
    //!
    //! \brief Clamps progress to [0,1] and returns the (rounded) fill
    //! position and percentage together, so every progress-bar variant
    //! agrees on where the bar ends and what percentage it reports.
    //!
    inline void progress_bar_metrics(
      double const progress, int const width, int & pos, int & pct
    )
    {
      double const p{ std::clamp( progress, 0.0, 1.0 ) };
      // std::round (not a raw static_cast<int>) avoids floating point
      // truncation artifacts, e.g. 0.7*50 evaluating to 34.999999... .
      pos = static_cast<int>( std::round( width * p ) );
      pct = static_cast<int>( std::round( 100.0 * p ) );
    }
  }  // namespace detail

  //!
  //! \brief Generates a text-based progress bar as a string.
  //!
  //! This function creates a progress bar in the form of a string
  //! representation, which visually indicates the progress of a task. The
  //! progress is shown using equal signs, a greater-than symbol for the current
  //! position, and underscores for the remaining part of the bar.
  //!
  //! \param progress A value between 0.0 and 1.0 representing the progress of
  //! the task (values outside this range are clamped).
  //! \param width    The total width of the progress bar (number of
  //! characters).
  //!
  //! \return A string representing the progress bar along with the percentage
  //! completed.
  //!
  inline string progress_bar( double const progress, int const width )
  {
    int pos, pct;
    detail::progress_bar_metrics( progress, width, pos, pct );
    string res{ "[" };
    for ( int i{ 0 }; i < width; ++i )
    {
      if ( i < pos )
        res += '=';
      else if ( i == pos )
        res += '>';
      else
        res += '_';
    }
    res += fmt::format( "] {:3d}%", pct );
    return res;
  }

  //!
  //! \brief Outputs an enhanced text-based progress bar to the specified output
  //! stream.
  //!
  //! This function generates a more detailed progress bar using various
  //! characters to represent different levels of completion. It outputs the
  //! progress bar directly to the provided output stream, including a
  //! percentage and an optional message.
  //!
  //! \param progress A value between 0.0 and 1.0 representing the progress of
  //! the task (values outside this range are clamped).
  //! \param width    The total width of the progress bar (number of
  //! characters).
  //! \param msg      An optional message to display alongside the progress bar.
  //!
  inline string progress_bar2( double const progress, int const width, string_view const msg )
  {
    static std::vector<string> const ch{ ".", "▏", "▎", "▍", "▌", "▋", "▊", "▉", "█" };
    double const p{ std::clamp( progress, 0.0, 1.0 ) };
    double const              ww{ width * p };
    int const                 pos{ static_cast<int>( floor( ww ) ) };
    int const                 frac8{ static_cast<int>( round( 8 * ( ww - pos ) ) ) };
    std::string               res = "[";
    for ( int i{ 0 }; i < width; ++i )
    {
      if ( i < pos )
        res += ch.back();
      else if ( i > pos )
        res += ch.front();
      else
        res += ch[frac8];
    }
    res += "] ";
    res += fmt::format( "{:3.0f}% {}", ceil( 100 * p ), msg );
    return res;
  }

  //!
  //! \brief Outputs a text-based progress bar to the specified output stream.
  //!
  //! Renders the exact same bar as \ref progress_bar (string overload) and
  //! writes it in place (using \\r plus an ANSI clear-to-end-of-line, so a
  //! shorter \c msg never leaves stray characters from a previous, longer
  //! call).
  //!
  //! \param s        The output stream where the progress bar will be
  //! displayed.
  //! \param progress A value between 0.0 and 1.0 representing the progress of
  //! the task (values outside this range are clamped).
  //! \param width    The total width of the progress bar (number of
  //! characters).
  //! \param msg      An optional message to display alongside the progress bar.
  //!
  inline void progress_bar( ostream_type & s, double const progress, int const width, string_view const msg )
  {
    s << progress_bar( progress, width ) << ' ' << msg << "\033[K\r" << std::flush;
  }

  //!
  //! \brief Outputs an enhanced text-based progress bar to the specified output
  //! stream.
  //!
  //! Renders the exact same bar as \ref progress_bar2 (string overload) and
  //! writes it in place (using \\r plus an ANSI clear-to-end-of-line, so a
  //! shorter \c msg never leaves stray characters from a previous, longer
  //! call).
  //!
  //! \param s        The output stream where the progress bar will be
  //! displayed.
  //! \param progress A value between 0.0 and 1.0 representing the progress of
  //! the task (values outside this range are clamped).
  //! \param width    The total width of the progress bar (number of
  //! characters).
  //! \param msg      An optional message to display alongside the progress bar.
  //!
  inline void progress_bar2( ostream_type & s, double const progress, int const width, string_view const msg )
  {
    s << progress_bar2( progress, width, msg ) << "\033[K\r" << std::flush;
  }

}  // namespace Utils

//
// eof: Utils_progress_bar.hh
//

#endif
