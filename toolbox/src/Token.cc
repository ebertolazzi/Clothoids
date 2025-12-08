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
// file: Token.cc
//

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include "Utils.hh"

namespace Utils
{
  //!
  //! \brief Retrieves the next token from the input string.
  //!
  //! This method scans the input string, starting from the current offset, and
  //! identifies the next substring that is not a delimiter. It updates the
  //! offset for the next call to ensure that the subsequent token is correctly
  //! located.
  //!
  //! A token is defined as a sequence of characters that is not contained in
  //! the specified set of delimiters. If no tokens are found, the method will
  //! return false. Otherwise, it will extract the token and return true.
  //!
  //! \return true if a token was found; false if there are no more tokens in
  //! the string.
  //!
  //! \note After calling this method, the extracted token can be accessed via
  //!       the member variable `m_token`.
  //!
  bool
  Tokenizer::next_token()
  {
    size_t const i{ m_string.find_first_not_of( m_delimiters, m_offset ) };
    m_offset = m_string.length();
    if ( string::npos == i ) return false;

    size_t const j{ m_string.find_first_of( m_delimiters, i ) };
    if ( string::npos == j )
    {
      m_token = m_string.substr( i );
      return true;
    }

    m_token  = m_string.substr( i, j - i );
    m_offset = j;
    return true;
  }

  //!
  //! \brief Splits a string into tokens based on specified delimiters.
  //!
  //! This function takes an input string and separates it into substrings
  //! (tokens) based on the specified delimiter string. The tokens are stored in
  //! the provided vector.
  //!
  //! \param str The input string to be split into tokens.
  //! \param sep A string containing the delimiter characters.
  //! \param arr A reference to a vector where the resulting tokens will be
  //! stored.
  //!
  //! This function uses the `Tokenizer` class to parse the input string and
  //! retrieve tokens until there are no more tokens left.
  //!
  //! \note The resulting tokens will be stored in the order they appear in the
  //! input string.
  //!
  void
  split_string( string const & str, string const & sep, vector<string> & arr )
  {
    Tokenizer tk( str, sep );
    while ( tk.next_token() ) arr.emplace_back( tk.get_token() );
  }

}  // namespace Utils

#endif

//
// eof: Token.cc
//
