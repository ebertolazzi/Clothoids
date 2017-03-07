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
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

//
// file: GenericContainerSupport.cc
//

#include "GenericContainer.hh"
#include <iomanip>

namespace GenericContainerNamespace {
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
  tokenizeString( std::string const & str,
                  vec_string_type   & tokens,
                  std::string const & delimiters ) {

    tokens.clear() ;

    // Skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while ( std::string::npos != pos || std::string::npos != lastPos ) {
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(lastPos, pos - lastPos));
      // Skip delimiters.  Note the "not_of"
      lastPos = str.find_first_not_of(delimiters, pos);
      // Find next "non-delimiter"
      pos = str.find_first_of(delimiters, lastPos);
    }
    // work around for end line delimiters
    if ( tokens.size() > 0 )
      if ( tokens.back()[0] == '\n' || tokens.back()[0] == '\r' )
        tokens.pop_back() ;
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
  getLineAndSkipComments( std::istream      & stream,
                          std::string       & line,
                          std::string const & commentchars ) {
    unsigned nl = 0 ;
    do {
      if ( stream.fail() ) return 0 ;
      getline( stream, line ) ;
      ++nl ;
    } while ( line.find_first_of(commentchars) != std::string::npos ) ; // ????????????
    return nl ;
  }

  // -------------------------------------------------------
  // original code by Francesco Biral (francesco.biral@unitn.it)
  GenericContainer const &
  GenericContainer::writeFormattedData( std::basic_ostream<char> & stream,
                                        char const delimiter ) const {
    GC_ASSERT( exists("headers"),
               "writeFormattedData, missing field `headers` in container") ;
    GC_ASSERT( exists("data"),
               "writeFormattedData, missing field `data` in container") ;
    GenericContainer const & data    = (*this)("data") ;
    vec_string_type  const & headers = (*this)("headers").get_vec_string(" writeFormattedData, `header` field must be `vec_string_type`") ;
    if ( (*this)("data").get_type() == GC_MAT_REAL )
      writeTable( headers, data.get_mat_real(), stream, delimiter ) ;
    else
      writeTable( headers, data.get_vector(), stream, delimiter ) ;
    return *this ;
  }

  // -------------------------------------------------------
  // original code by Francesco Biral (francesco.biral@unitn.it)
  GenericContainer &
  GenericContainer::readFormattedData( std::istream & stream,
                                       char const commentChars[],
                                       char const delimiters[] ) {
    //read a line
    std::string line ;

    this -> set_map() ;
    GenericContainer & tmp = (*this)["headers"] ;
    
    std::cout << tmp.get_type_name() << '\n';
    vec_string_type & headers = tmp.set_vec_string() ;

    // reading header line
    unsigned nline = getLineAndSkipComments( stream, line, commentChars ) ; // read  line
    tokenizeString( line, headers, delimiters ) ; // tokenize line
    unsigned ncol = unsigned( headers.size() ) ;

    vector_type & data = (*this)["data"].set_vector(ncol) ;
    for ( unsigned icol = 0 ; icol < ncol ; ++icol ) data[icol].set_vec_real() ;

    // read data by line
    unsigned nread ;
    vec_string_type tokens ;
    while ( (nread=getLineAndSkipComments( stream, line, commentChars )) > 0 ) {
      nline += nread ;

      // read line and convert into vector of strings
      tokenizeString( line, tokens, delimiters ) ;
      if ( tokens.size() == 0 ) break ; // riga vuota!

      GC_ASSERT( unsigned(tokens.size()) == ncol,
                 "readFormattedDataFile, in reading line: " << nline <<
                 " expected " << ncol << " found: " << tokens.size() ) ;

      // store data in row vector
      for ( unsigned icol = 0 ; icol < ncol ; ++icol )
        data[icol].get_vec_real().push_back(atof(tokens[icol].c_str()) );
    }
    return *this ;
  }
}

//
// eof: GenericContainerSupport.cc
//
