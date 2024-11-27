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

//
// https://stackoverflow.com/questions/53849/how-do-i-tokenize-a-string-in-c
//

namespace Utils {

  using std::size_t;
  using std::string;
  using std::vector;

  class Tokenizer {
  protected:
    size_t            m_offset;
    std::string const m_string;
    std::string       m_token;
    std::string const m_delimiters;
  public:
    Tokenizer(
      std::string str,
      std::string delimiters
    )
    : m_offset(0)
    , m_string(std::move(str))
    , m_delimiters(std::move(delimiters))
    { }

    std::string get_token() const { return m_token; }

    bool next_token();
  };

  void
  split_string(
    std::string const        & str,
    std::string const        & sep,
    std::vector<std::string> & arr
  );

}

///
/// eof: Token.hxx
///
