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

///
/// file: Token.hxx
///

#pragma once

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#ifndef UTILS_TOKEN_dot_HXX
#define UTILS_TOKEN_dot_HXX
#endif

//
// https://stackoverflow.com/questions/53849/how-do-i-tokenize-a-string-in-c
//

namespace Utils {

  class Tokenizer {
  protected:
    size_t       m_offset;
    string const m_string;
    string       m_token;
    string const m_delimiters;
  public:
    Tokenizer( string const & str, string const & delimiters )
    : m_offset(0)
    , m_string(str)
    , m_delimiters(delimiters)
    { }

    string get_token() const { return m_token; }

    bool next_token();
  };

  void
  split_string(
    string const   & str,
    string const   & sep,
    vector<string> & arr
  );

}


#endif

///
/// eof: Token.hxx
///
