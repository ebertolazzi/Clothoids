/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2023                                                      |
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

///
/// file: Utils_progress_bar.cc
///

#include "Utils.hh"

namespace Utils {

  using std::ceil;
  using std::floor;
  using std::round;

  string
  progress_bar( double progress, int width ) {
    string res{"["};
    int pos = int(width * progress);
    for (int i = 0; i < width; ++i) {
      if      (i  < pos) res += '=';
      else if (i == pos) res += '>';
      else               res += '_';
    }
    res += fmt::format( "] {:3.0f}%", ceil(100*progress) );
    return res;
  }

  void
  progress_bar( ostream & s, double progress, int width, char const * msg ) {
    s << "[";
    int pos = int(width * progress);
    for (int i = 0; i < width; ++i)
      s << (i<pos? u8"\u25A0":u8"\u25A1");
    fmt::print( s, "] {:3.0f}% {}\r", ceil(100*progress), msg );
    s << std::flush;
  }

  void
  progress_bar2( ostream & s, double progress, int width, char const * msg ) {
    vector<string> ch{"_", "▏", "▎", "▍", "▌", "▋", "▊", "▉", "█"};
    s << "[";
    double ww = width * progress;
    int pos   = int(floor(ww));
    int frac8 = int(round(8*(ww-pos)));
    for (int i = 0; i < width; ++i) {
      if      (i < pos) s << ch.back();
      else if (i > pos) s << ' ';
      else              s << ch[frac8];
    }
    fmt::print( s, "] {:3.0f}% {}\r", ceil(100*progress), msg );
    s << std::flush;
  }

}

///
/// eof: Utils_progress_bar.c
///
