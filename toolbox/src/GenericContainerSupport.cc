/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2013                                                      |
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
// file: GenericContainerSupport.cc
//

#include "GenericContainer/GenericContainer.hh"
#include <iomanip>
#include <fstream>
#include <cstdlib>

namespace GC_namespace {

  using std::vector;
  using std::string;
  using std::atof;
  using std::strtod;

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  /*
  //   _____     _              _
  //  |_   _|__ | | _____ _ __ (_)_______
  //    | |/ _ \| |/ / _ \ '_ \| |_  / _ \
  //    | | (_) |   <  __/ | | | |/ /  __/
  //    |_|\___/|_|\_\___|_| |_|_/___\___|
  */
  // Tokenizes a string parsing the string using space delimiter (default) or given one.
  // The string is slitted into substrings stored into a vector of strings.
  //
  // Original code by Francesco Biral (francesco.biral@unitn.it)
  static
  void
  tokenizeString(
    string_view const str,
    vec_string_type & tokens,
    string_view const delimiters
  ) {

    tokens.clear();

    // Skip delimiters at beginning.
    string_type::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string_type::size_type pos = str.find_first_of(delimiters, lastPos);

    while ( string_type::npos != pos || string_type::npos != lastPos ) {
      // Found a token, add it to the vector.
      tokens.emplace_back(str.substr(lastPos, pos - lastPos));
      // Skip delimiters.  Note the "not_of"
      lastPos = str.find_first_not_of(delimiters, pos);
      // Find next "non-delimiter"
      pos = str.find_first_of(delimiters, lastPos);
    }
    // work around for end line delimiters
    if ( !tokens.empty() )
      if ( tokens.back()[0] == '\n' || tokens.back()[0] == '\r' )
        tokens.pop_back();
  }

  /*
  //   ___    _____          _   _ _
  //  |_ _|  / / _ \   _   _| |_(_) |___
  //   | |  / / | | | | | | | __| | / __|
  //   | | / /| |_| | | |_| | |_| | \__ \
  //  |___/_/  \___/   \__,_|\__|_|_|___/
  */

  static
  unsigned
  get_line_and_skip_comments(
    istream_type    & stream,
    string_type     & line,
    string_view const commentchars
  ) {
    unsigned nl = 0;
    do {
      if ( stream.fail() ) return 0;
      getline( stream, line );
      ++nl;
    } while ( line.find_first_of(commentchars) != string_type::npos ); // ????????????
    return nl;
  }

  static
  unsigned
  get_line_and_skip_comments2(
    istream_type     & stream,
    string_type      & line,
    string_view const  commentchars,
    GenericContainer * ptr_pars
  ) {
    unsigned nl      = 0;
    bool     comment = true;
    while ( comment ) {
      if ( stream.fail() ) return 0;
      getline( stream, line );
      comment = line.empty() ||
                line.find_first_of(commentchars) != string_type::npos;
      if ( ptr_pars != nullptr && comment && line.length() > 2 && line[1] == '!' ) {
        vec_string_type a_eq_b;
        tokenizeString( line.substr(2), a_eq_b, "=" );
        if ( a_eq_b.size() > 1 ) {
          // parse parameter value = data
          // trim trailing spaces
          string lhs = a_eq_b[0];
          string rhs = a_eq_b[1];
          // trim
          if ( size_t const endpos{ lhs.find_last_not_of(" \t") };    string::npos != endpos   ) lhs = lhs.substr( 0, endpos+1 );
          if ( size_t const startpos{ lhs.find_first_not_of(" \t") }; string::npos != startpos ) lhs = lhs.substr( startpos );
          // trim
          if ( size_t const endpos{ rhs.find_last_not_of(" \t") };    string::npos != endpos   ) rhs = rhs.substr( 0, endpos+1 );
          if ( size_t const startpos{ rhs.find_first_not_of(" \t") }; string::npos != startpos ) rhs = rhs.substr( startpos );
          char * ptr;
          (*ptr_pars)[lhs] = strtod(rhs.data(),&ptr);
          if ( ptr == rhs.data() ) (*ptr_pars)[lhs] = rhs;
        }
      }
      ++nl;
    }
    return nl;
  }

  #endif

  // -------------------------------------------------------
  // original code by Francesco Biral (francesco.biral@unitn.it)
  GenericContainer const &
  GenericContainer::write_formatted_data(
    ostream_type & stream,
    char const     delimiter
  ) const {
    GC_ASSERT(
      this->exists("headers"),
      "write_formatted_data, missing field `headers` in container"
    );
    GC_ASSERT(
      this->exists("data"),
      "write_formatted_data, missing field `data` in container"
    );
    GenericContainer const & data    = (*this)("data");
    vec_string_type  const & headers = (*this)("headers").get_vec_string("write_formatted_data, `header` field must be `vec_string_type`");
    if ( (*this)("data").get_type() == GC_type::MAT_REAL )
      write_table( headers, data.get_mat_real(), stream, delimiter );
    else
      write_table( headers, data.get_vector(), stream, delimiter );
    return *this;
  }

  // -------------------------------------------------------
  // original code by Francesco Biral (francesco.biral@unitn.it)
  GenericContainer &
  GenericContainer::read_formatted_data(
    istream_type & stream,
    char const     commentChars[],
    char const     delimiters[]
  ) {
    //read a line
    string_type line;

    this->set_map();
    GenericContainer & tmp = (*this)["headers"];
    vec_string_type & headers = tmp.set_vec_string();

    // reading header line
    unsigned nline { get_line_and_skip_comments( stream, line, commentChars ) }; // read  line
    tokenizeString( line, headers, delimiters ); // tokenize line
    unsigned const ncol { static_cast<unsigned>(headers.size()) };

    vector_type & data { (*this)["data"].set_vector(ncol) };
    for ( unsigned icol{0}; icol < ncol; ++icol ) data[icol].set_vec_real();

    // read data by line
    unsigned nread;
    vec_string_type tokens;
    while ( (nread=get_line_and_skip_comments( stream, line, commentChars )) > 0 ) {
      nline += nread;

      // read line and convert into vector of strings
      tokenizeString( line, tokens, delimiters );
      if ( tokens.empty() ) break; // riga vuota!

      GC_ASSERT(
        static_cast<unsigned>(tokens.size()) == ncol,
        "read_formatted_data, in reading line: " << nline <<
        " expected " << ncol << " found: " << tokens.size()
      );

      // store data in row vector
      for ( unsigned icol = 0; icol < ncol; ++icol )
        data[icol].get_vec_real().push_back(atof(tokens[icol].data()) );
    }
    return *this;
  }

  // -------------------------------------------------------
  // original code by Francesco Biral (francesco.biral@unitn.it)
  GenericContainer &
  GenericContainer::read_formatted_data2(
    istream_type   & stream,
    char const       commentChars[],
    char const       delimiters[],
    GenericContainer ptr_pars[]
  ) {
    //read a line
    string_type line;

    this->set_map();
    GenericContainer & tmp = (*this)["headers"];
    vec_string_type & headers = tmp.set_vec_string();

    // reading header line
    auto nline{ get_line_and_skip_comments2( stream, line, commentChars, ptr_pars ) }; // read  line
    tokenizeString( line, headers, delimiters ); // tokenize line
    auto const ncol { static_cast<unsigned>(headers.size()) };

    vector<vec_real_type*> pcolumns(ncol);

    GenericContainer & data = (*this)["data"];
    for ( unsigned icol = 0; icol < ncol; ++icol ) {
      GenericContainer & ICOL = data[headers[icol]];
      ICOL.set_vec_real();
      pcolumns[icol] = &ICOL.get_vec_real();
    }

    // read data by line
    unsigned nread;
    vec_string_type tokens;
    while ( (nread=get_line_and_skip_comments2( stream, line, commentChars, ptr_pars )) > 0 ) {
      nline += nread;

      // read line and convert into vector of strings
      tokenizeString( line, tokens, delimiters );
      if ( tokens.empty() ) break; // riga vuota!

      GC_ASSERT(
        static_cast<unsigned>(tokens.size()) == ncol,
        "read_formatted_data2, in reading line: " << nline <<
        " expected " << ncol << " found: " << tokens.size()
      );

      // store data in row vector
      for ( unsigned icol = 0; icol < ncol; ++icol )
        pcolumns[icol]->push_back(atof(tokens[icol].data()) );
    }
    return *this;
  }

  GenericContainer &
  GenericContainer::read_formatted_data(
    char const fname[],
    char const commentChars[],
    char const delimiters[]
  ) {
    std::ifstream file( fname );
    GC_ASSERT(
      file.good(),
      "read_formatted_data, failed to open file: ``" << fname << "''"
    )
    return read_formatted_data( file, commentChars, delimiters );
  }

  GenericContainer &
  GenericContainer::read_formatted_data2(
    char const       fname[],
    char const       commentChars[],
    char const       delimiters[],
    GenericContainer ptr_pars[]
  ) {
    std::ifstream file( fname );
    GC_ASSERT(
      file.good(),
      "read_formatted_data2, failed to open file: ``" << fname << "''"
    )
    return read_formatted_data2( file, commentChars, delimiters, ptr_pars );
  }

}

//
// eof: GenericContainerSupport.cc
//
