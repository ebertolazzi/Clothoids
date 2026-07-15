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

#ifdef _MSC_VER
#pragma warning( disable : 4661 )
#endif

#include "GenericContainer/GenericContainer.hh"
#include <iomanip>
#include <fstream>
#include <cstdlib>

namespace GC_namespace
{

  using std::atof;
  using std::string;
  using std::strtod;
  using std::vector;


  bool GenericContainer::from_file( string_view file_name )
  {
    string_type const path{ file_name };
    std::ifstream     file( path );
    if ( !file.is_open() ) return false;

    // Utility to check file extension
    auto ends_with = []( string_view str, string_view suffix ) -> bool
    { return str.size() >= suffix.size() && str.compare( str.size() - suffix.size(), suffix.size(), suffix ) == 0; };

    // JSON support was dropped from the core library: nlohmann::json (via
    // include/GenericContainer/GenericContainerInterface_nlohmann.hh) is the
    // replacement for C++ callers -- j.get<GenericContainer>() / j = gc --
    // but there's no FFI/file-dispatch equivalent, so a ".json" file is
    // simply not handled here.
    //
    // from_yaml/from_toml ARE implemented, but as separate opt-in libraries
    // (src_yaml_interface/, src_toml_interface/, built as GenericContainer::
    // Yaml / GenericContainer::Toml -- see their CMakeLists.txt) rather than
    // linked into this core library, precisely so the core stays a
    // self-contained, standalone-linkable target (a shared build must
    // resolve every symbol its object files reference at link time; calling
    // through to a backend this library doesn't itself link would break
    // that regardless of whether the branch is ever taken). Callers who
    // need yaml/toml should link the relevant optional library and call
    // gc.from_yaml()/gc.from_toml() directly instead of going through this
    // extension-sniffing convenience dispatch.
    if ( ends_with( file_name, ".yaml" ) || ends_with( file_name, ".yml" ) || ends_with( file_name, ".toml" ) )
    {
      file.close();
      GC_assert(
        false,
        "GenericContainer::from_file: yaml/toml support is not linked into this library; "
        "link GenericContainer::Yaml or GenericContainer::Toml and call from_yaml()/from_toml() directly" );
    }
    file.close();
    return false;
  }

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
  static void tokenizeString( string_view const str, vec_string_type & tokens, string_view const delimiters )
  {
    tokens.clear();

    // Skip delimiters at beginning.
    string_type::size_type lastPos = str.find_first_not_of( delimiters, 0 );
    // Find first "non-delimiter".
    string_type::size_type pos = str.find_first_of( delimiters, lastPos );

    while ( string_type::npos != pos || string_type::npos != lastPos )
    {
      // Found a token, add it to the vector.
      tokens.emplace_back( str.substr( lastPos, pos - lastPos ) );
      // Skip delimiters.  Note the "not_of"
      lastPos = str.find_first_not_of( delimiters, pos );
      // Find next "non-delimiter"
      pos = str.find_first_of( delimiters, lastPos );
    }
    // work around for end line delimiters
    if ( !tokens.empty() )
      if ( tokens.back()[0] == '\n' || tokens.back()[0] == '\r' ) tokens.pop_back();
  }

  /*
  //   ___    _____          _   _ _
  //  |_ _|  / / _ \   _   _| |_(_) |___
  //   | |  / / | | | | | | | __| | / __|
  //   | | / /| |_| | | |_| | |_| | \__ \
  //  |___/_/  \___/   \__,_|\__|_|_|___/
  */

  static std::size_t get_line_and_skip_comments(
    istream_type &    stream,
    string_type &     line,
    string_view const commentchars )
  {
    std::size_t nl = 0;
    do
    {
      if ( stream.fail() ) return 0;
      getline( stream, line );
      ++nl;
    } while ( line.find_first_of( commentchars ) != string_type::npos );  // ????????????
    return nl;
  }

  static std::size_t get_line_and_skip_comments2(
    istream_type &     stream,
    string_type &      line,
    string_view const  commentchars,
    GenericContainer * ptr_pars )
  {
    std::size_t nl      = 0;
    bool        comment = true;
    while ( comment )
    {
      if ( stream.fail() ) return 0;
      getline( stream, line );
      comment = line.empty() || line.find_first_of( commentchars ) != string_type::npos;
      if ( ptr_pars != nullptr && comment && line.length() > 2 && line[1] == '!' )
      {
        vec_string_type a_eq_b;
        tokenizeString( line.substr( 2 ), a_eq_b, "=" );
        if ( a_eq_b.size() > 1 )
        {
          // parse parameter value = data
          // trim trailing spaces
          string lhs = a_eq_b[0];
          string rhs = a_eq_b[1];
          // trim
          if ( size_t const endpos{ lhs.find_last_not_of( " \t" ) }; string::npos != endpos )
            lhs = lhs.substr( 0, endpos + 1 );
          if ( size_t const startpos{ lhs.find_first_not_of( " \t" ) }; string::npos != startpos )
            lhs = lhs.substr( startpos );
          // trim
          if ( size_t const endpos{ rhs.find_last_not_of( " \t" ) }; string::npos != endpos )
            rhs = rhs.substr( 0, endpos + 1 );
          if ( size_t const startpos{ rhs.find_first_not_of( " \t" ) }; string::npos != startpos )
            rhs = rhs.substr( startpos );
          char * ptr;
          ( *ptr_pars )[lhs] = strtod( rhs.data(), &ptr );
          if ( ptr == rhs.data() ) ( *ptr_pars )[lhs] = rhs;
        }
      }
      ++nl;
    }
    return nl;
  }

#endif

  // -------------------------------------------------------
  // original code by Francesco Biral (francesco.biral@unitn.it)
  GenericContainer const & GenericContainer::write_formatted_data( ostream_type & stream, char const delimiter ) const
  {
    GC_assert( this->exists( "headers" ), "write_formatted_data, missing field `headers` in container" );
    GC_assert( this->exists( "data" ), "write_formatted_data, missing field `data` in container" );
    GenericContainer const & data = ( *this )( "data" );
    vec_string_type const &  headers =
      ( *this )( "headers" ).get_vec_string( "write_formatted_data, `header` field must be `vec_string_type`" );
    if ( ( *this )( "data" ).get_type() == GC_type::MAT_REAL )
      write_table( headers, data.get_mat_real(), stream, delimiter );
    else
      write_table( headers, data.get_vector(), stream, delimiter );
    return *this;
  }

  // -------------------------------------------------------
  // original code by Francesco Biral (francesco.biral@unitn.it)
  GenericContainer & GenericContainer::read_formatted_data(
    istream_type & stream,
    char const     commentChars[],
    char const     delimiters[] )
  {
    // read a line
    string_type line;

    this->set_map();
    GenericContainer & tmp     = ( *this )["headers"];
    vec_string_type &  headers = tmp.set_vec_string();

    // reading header line
    std::size_t nline{ get_line_and_skip_comments( stream, line, commentChars ) };  // read  line
    tokenizeString( line, headers, delimiters );                                    // tokenize line
    std::size_t const ncol{ static_cast<std::size_t>( headers.size() ) };

    vector_type & data{ ( *this )["data"].set_vector( ncol ) };
    for ( std::size_t icol{ 0 }; icol < ncol; ++icol ) data[icol].set_vec_real();

    // read data by line
    std::size_t     nread;
    vec_string_type tokens;
    while ( ( nread = get_line_and_skip_comments( stream, line, commentChars ) ) > 0 )
    {
      nline += nread;

      // read line and convert into vector of strings
      tokenizeString( line, tokens, delimiters );
      if ( tokens.empty() ) break;  // riga vuota!

      GC_assert(
        static_cast<std::size_t>( tokens.size() ) == ncol,
        "read_formatted_data, in reading line: {} expected {} found: {}",
        nline,
        ncol,
        tokens.size() );

      // store data in row vector
      for ( std::size_t icol = 0; icol < ncol; ++icol )
        data[icol].get_vec_real().push_back( atof( tokens[icol].data() ) );
    }
    return *this;
  }

  // -------------------------------------------------------
  // original code by Francesco Biral (francesco.biral@unitn.it)
  GenericContainer & GenericContainer::read_formatted_data2(
    istream_type &   stream,
    char const       commentChars[],
    char const       delimiters[],
    GenericContainer ptr_pars[] )
  {
    // read a line
    string_type line;

    this->set_map();
    GenericContainer & tmp     = ( *this )["headers"];
    vec_string_type &  headers = tmp.set_vec_string();

    // reading header line
    auto nline{ get_line_and_skip_comments2( stream, line, commentChars, ptr_pars ) };  // read  line
    tokenizeString( line, headers, delimiters );                                        // tokenize line
    auto const ncol{ static_cast<std::size_t>( headers.size() ) };

    vector<vec_real_type *> pcolumns( ncol );

    GenericContainer & data = ( *this )["data"];
    for ( std::size_t icol = 0; icol < ncol; ++icol )
    {
      GenericContainer & ICOL = data[headers[icol]];
      ICOL.set_vec_real();
      pcolumns[icol] = &ICOL.get_vec_real();
    }

    // read data by line
    std::size_t     nread;
    vec_string_type tokens;
    while ( ( nread = get_line_and_skip_comments2( stream, line, commentChars, ptr_pars ) ) > 0 )
    {
      nline += nread;

      // read line and convert into vector of strings
      tokenizeString( line, tokens, delimiters );
      if ( tokens.empty() ) break;  // riga vuota!

      GC_assert(
        static_cast<std::size_t>( tokens.size() ) == ncol,
        "read_formatted_data2, in reading line: {} expected {} found: {}",
        nline,
        ncol,
        tokens.size() );

      // store data in row vector
      for ( std::size_t icol = 0; icol < ncol; ++icol ) pcolumns[icol]->push_back( atof( tokens[icol].data() ) );
    }
    return *this;
  }

  GenericContainer & GenericContainer::read_formatted_data(
    char const fname[],
    char const commentChars[],
    char const delimiters[] )
  {
    std::ifstream file( fname );
    GC_assert( file.good(), "read_formatted_data, failed to open file: ``{}''", fname );
    return read_formatted_data( file, commentChars, delimiters );
  }

  GenericContainer & GenericContainer::read_formatted_data2(
    char const       fname[],
    char const       commentChars[],
    char const       delimiters[],
    GenericContainer ptr_pars[] )
  {
    std::ifstream file( fname );
    GC_assert( file.good(), "read_formatted_data2, failed to open file: ``{}''", fname );
    return read_formatted_data2( file, commentChars, delimiters, ptr_pars );
  }

}  // namespace GC_namespace

//
// eof: GenericContainerSupport.cc
//
