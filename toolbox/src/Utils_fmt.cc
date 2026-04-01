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
// file: Utils_fmt.cc
//

#include "Utils_fmt.hh"

namespace Utils
{

  string fmt_table_row(
    unsigned const    width,
    string_view const L,
    string_view const R,
    string_view const F,
    string_view const title,
    string_view const align )
  {
    string const FMT{ fmt::format( "{{}}{{:{}{}{}}}{{}}", F, align, width ) };
    return fmt::format( fmt::runtime( FMT ), L, title, R );
  }

  string fmt_table_row(
    unsigned const                           width,
    string_view const                        L,
    string_view const                        C,
    string_view const                        R,
    string_view const                        F,
    std::initializer_list<string_view> const names,
    string_view                              align )
  {
    unsigned const N{ static_cast<unsigned>( names.size() ) };
    unsigned const ww{ width + 1 - N };
    unsigned const w{ ww / N };
    unsigned const r{ ww - w * N };
    unsigned const r2{ r / 2 };
    UTILS_ASSERT( w > 3, "fmt_table_row( width={}, ... ) no space to print\n", width );
    string const FMT{ fmt::format( "{}{{:{}{}{}}}{}", F, F, align, w - 2, F ) };
    string       res{ L };
    unsigned     k{ 0 };
    while ( k < r2 )
    {
      res += ' ';
      ++k;
    }
    unsigned i{ 0 };
    for ( auto n : names )
    {
      res += fmt::format( fmt::runtime( FMT ), n );
      if ( ++i != N ) res += C;
    }
    while ( k < r )
    {
      res += ' ';
      ++k;
    }
    res += R;
    return res;
  }

  string fmt_table_row(
    unsigned const    width,
    string_view const L,
    string_view const C,
    string_view const R,
    string_view const F,
    unsigned const    N )
  {
    unsigned const ww{ width + 1 - N };
    unsigned const w{ ww / N };
    unsigned const r{ ww - w * N };
    unsigned const r2{ r / 2 };
    UTILS_ASSERT( w > 3, "fmt_table_row( width={}, ... ) no space to print\n", width );

    string str;
    for ( unsigned i{ 0 }; i < w; ++i ) str += F;

    string   res{ L };
    unsigned k{ 0 };
    while ( k < r2 )
    {
      res += F;
      ++k;
    }
    for ( unsigned i{ 1 }; i <= N; ++i )
    {
      res += str;
      if ( i != N ) res += C;
    }
    while ( k < r )
    {
      res += F;
      ++k;
    }
    res += R;
    return res;
  }

}  // namespace Utils

//
// EOF: Utils_fmt.cc
//
