/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2022                                                      |
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
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

//
// file: Utils_fmt.cc
//

#include "Utils_fmt.hh"

namespace Utils {

  string
  fmt_table_row(
    unsigned    width,
    string_view L,
    string_view R,
    string_view F,
    string_view title,
    string_view align
  ) {
    string FMT{ fmt::format( "{{}}{{:{}{}{}}}{{}}", F, align, width ) };
    return fmt::format( FMT, L, title, R );
  }
  
  string
  fmt_table_row(
    unsigned    width,
    string_view L,
    string_view C,
    string_view R,
    string_view F,
    std::initializer_list<string_view> names,
    string_view align
  ) {
    unsigned N  { unsigned(names.size()) };
    unsigned ww { width+1-N };
    unsigned w  { ww/N };
    unsigned r  { ww-w*N };
    unsigned r2 { r/2 };
    UTILS_ASSERT( w > 3, "fmt_table_row( width={}, ... ) no space to print\n", width );
    string FMT{ fmt::format("{}{{:{}{}{}}}{}", F, F, align, w-2, F ) };
    string res{ L };
    unsigned k{0};
    while ( k < r2 ) { res += ' '; ++k; }
    unsigned i{0};
    for ( auto n : names ) {
      res += fmt::format( FMT, n );
      if ( ++i != N ) res += C;
    };
    while ( k < r ) { res += ' '; ++k; }
    res += R;
    return res;
  }
  
  string
  fmt_table_row(
    unsigned    width,
    string_view L,
    string_view C,
    string_view R,
    string_view F,
    unsigned    N
  ) {
    unsigned ww { width+1-N };
    unsigned w  { ww/N };
    unsigned r  { ww-w*N };
    unsigned r2 { r/2 };
    UTILS_ASSERT( w > 3, "fmt_table_row( width={}, ... ) no space to print\n", width );

    string str{""};
    for ( unsigned i{0}; i < w; ++i ) str += F;

    string res{ L };
    unsigned k{0};
    while ( k < r2 ) { res += F; ++k; }
    for ( unsigned i{1}; i <= N; ++i ) { res += str; if ( i != N ) res += C; };
    while ( k < r ) { res += F; ++k; }
    res += R;
    return res;
  }
  
}

//
// EOF: Utils_fmt.cc
//
