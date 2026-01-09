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
// file: Token.hxx (Header-Only Version)
//

#ifndef UTILS_TOKEN_HH
#define UTILS_TOKEN_HH

#include "Utils.hh"

namespace Utils
{
  using std::size_t;
  using std::string;
  using std::vector;

  /*!
   * \addtogroup UTILS
   * @{
   */

  //!!
  //! \brief A utility class for tokenizing strings.
  //!
  //! The `Tokenizer` class provides a simple and efficient mechanism for
  //! breaking a string into smaller components (tokens) based on specified
  //! delimiter(s). It facilitates sequential token extraction from the input
  //! string and supports custom delimiters.
  //!
  //! **Usage Example:**
  //! \code{cpp}
  //! Utils::Tokenizer tokenizer("Hello, world! Welcome to Tokenizer.", " ,!");
  //! while (tokenizer.next_token()) {
  //!   std::cout << tokenizer.get_token() << std::endl;
  //! }
  //! \endcode
  //!
  class Tokenizer
  {
  protected:
    size_t       m_offset;      //!< Current position in the string.
    string const m_string;      //!< The original string to tokenize.
    string       m_token;       //!< The current token extracted.
    string const m_delimiters;  //!< Characters used to separate tokens.
  public:
    //!
    //! Constructs a Tokenizer object with the provided string and delimiters.
    //!
    //! \param str The input string to tokenize.
    //! \param delimiters A string containing delimiter characters.
    //!
    Tokenizer( string str, string delimiters )
      : m_offset( 0 ), m_string( std::move( str ) ), m_delimiters( std::move( delimiters ) )
    {
    }

    //!
    //! Retrieves the last extracted token.
    //!
    //! \return The last token that was retrieved using `next_token()`.
    //!
    string const & get_token() const { return m_token; }

    //!
    //! Advances the tokenizer to the next token in the string.
    //!
    //! This function searches the string for the next token, updating the
    //! internal state of the `Tokenizer` object. If a token is found, it sets
    //! the `m_token` member variable to the new token and returns `true`. If no
    //! more tokens are available, it returns `false`.
    //!
    //! \return `true` if a new token was found; `false` if no more tokens are
    //! available.
    //!
    bool next_token()
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
  };

  //!
  //! Splits a string into tokens based on specified delimiters.
  //!
  //! This function takes an input string and splits it into a vector of tokens,
  //! using the specified separators to identify where to split the string.
  //!
  //! \param str The input string to split.
  //! \param sep A string containing delimiter characters.
  //! \param arr A vector that will be filled with the resulting tokens.
  //!
  //! \note The resulting tokens are stored in the provided vector, and the
  //! vector
  //!       will be resized to accommodate all tokens extracted from the input
  //!       string.
  //!
  inline void split_string( string const & str, string const & sep, vector<string> & arr )
  {
    Tokenizer tk( str, sep );
    while ( tk.next_token() ) arr.emplace_back( tk.get_token() );
  }

  /*! @} */

}  // namespace Utils

#endif  // UTILS_TOKEN_HXX
