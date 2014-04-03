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
// file: GenericContainerTables.cc
//

#include "GenericContainer.hh"
#include <iomanip>
#include <string>

namespace GC {

  using namespace std ;

  void
  writeTable( vec_string_type const    & headers,
              vector_type     const    & data,
              std::basic_ostream<char> & stream,
              char const                 delimiter ) {

    unsigned ncol = (unsigned)headers.size() ;
    unsigned nrow = (unsigned)data[0].get_num_elements() ;

    stream << headers[0] ;
    for ( unsigned icol = 1 ; icol < ncol ; ++icol )
      stream << delimiter << headers[icol] ;
    stream << '\n' ;

    for ( unsigned row = 0 ; row < nrow ; ++row ) {
      stream << data[0].get_number(row) ;
      for ( unsigned icol = 1 ; icol < ncol ; ++icol )
        stream << delimiter << data[icol].get_number(row)  ;
      stream << '\n' ;
    }
  }

  void
  writeTableFormatted( vec_string_type const    & headers,
                       vector_type     const    & data,
                       std::basic_ostream<char> & stream ) {

    unsigned ncol = unsigned(headers.size()) ;
    unsigned nrow = data[0].get_num_elements() ;

    if ( ncol == 0 ) return ;

    // calcolo lunghezza massima stringhe headers
    unsigned ml = 0 ;
    vec_string_type::const_iterator is = headers.begin() ;
    for ( ; is != headers.end() ; ++is )
      if ( ml < is->length() ) ml = unsigned(is->length()) ;
    // taglio a lunghezza min/max
    if      ( ml < 8  ) ml = 8 ;
    else if ( ml > 20 ) ml = 20 ;

    string line = std::string(ncol*(ml+1)-1, '-') ;

    is = headers.begin() ;
    stream << std::setw(int(ml)) << *is ;
    for ( ++is ; is != headers.end() ; ++is )
      stream << " " << std::setw(int(ml)) << *is ;
    stream << '\n' << line << '\n';

    for ( unsigned row = 0 ; row < nrow ; ++row ) {
      stream << std::setw(int(ml)) << data[0].get_number(row) ;
      for ( unsigned icol = 1 ; icol < ncol ; ++icol )
        stream << " " << std::setw(int(ml)) << data[icol].get_number(row)  ;
      stream << '\n' ;
    }
    stream << line << '\n';
  }


}

//
// eof: GenericContainerTables.cc
//
