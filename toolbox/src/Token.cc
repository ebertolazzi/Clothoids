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
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include "Utils.hh"

namespace Utils {

  bool
  Tokenizer::next_token() {
    size_t i = m_string.find_first_not_of( m_delimiters, m_offset );
    m_offset = m_string.length();
    if ( string::npos == i ) return false;

    size_t j = m_string.find_first_of( m_delimiters, i );
    if ( string::npos == j ) {
      m_token = m_string.substr(i);
      return true;
    }

    m_token  = m_string.substr(i, j - i);
    m_offset = j;
    return true;
  }

  void
  split_string(
    string const   & str,
    string const   & sep,
    vector<string> & arr
  ) {
    Tokenizer tk( str, sep );
    while ( tk.next_token() ) arr.push_back(tk.get_token());
  }

}

#endif

///
/// eof: Token.cc
///

