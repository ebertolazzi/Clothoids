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

//
// file: Utils_progress_bar.cc
//

#include "Utils.hh"
#include "Utils_fmt.hh"

namespace Utils {

  using std::ceil;
  using std::floor;
  using std::round;

  //!
  //! \brief Generates a text-based progress bar as a string.
  //!
  //! This function creates a progress bar in the form of a string representation,
  //! which visually indicates the progress of a task. The progress is shown using
  //! equal signs, a greater-than symbol for the current position, and underscores
  //! for the remaining part of the bar.
  //!
  //! \param progress A value between 0.0 and 1.0 representing the progress of the task.
  //! \param width    The total width of the progress bar (number of characters).
  //!
  //! \return A string representing the progress bar along with the percentage completed.
  //!
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

  //!
  //! \brief Outputs a text-based progress bar to the specified output stream.
  //!
  //! This function generates a progress bar and outputs it directly to the provided
  //! output stream, indicating the current progress of a task. The visual representation
  //! uses filled and unfilled squares to indicate progress and includes a percentage
  //! along with an optional message.
  //!
  //! \param s        The output stream where the progress bar will be displayed.
  //! \param progress A value between 0.0 and 1.0 representing the progress of the task.
  //! \param width    The total width of the progress bar (number of characters).
  //! \param msg      An optional message to display alongside the progress bar.
  //!
  void
  progress_bar( ostream & s, double progress, int width, char const * msg ) {
    s << "[";
    int pos = int(width * progress);
    for (int i = 0; i < width; ++i)
      s << (i<pos? u8"\u25A0":u8"\u25A1");
    fmt::print( s, "] {:3.0f}% {}\r", ceil(100*progress), msg );
    s << std::flush;
  }

  //!
  //! \brief Outputs an enhanced text-based progress bar to the specified output stream.
  //!
  //! This function generates a more detailed progress bar using various characters to
  //! represent different levels of completion. It outputs the progress bar directly to
  //! the provided output stream, including a percentage and an optional message.
  //!
  //! \param s        The output stream where the progress bar will be displayed.
  //! \param progress A value between 0.0 and 1.0 representing the progress of the task.
  //! \param width    The total width of the progress bar (number of characters).
  //! \param msg      An optional message to display alongside the progress bar.
  //!
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

//
// eof: Utils_progress_bar.c
//
