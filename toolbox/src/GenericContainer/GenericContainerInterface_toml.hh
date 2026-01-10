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

#pragma once

#ifndef GENERIC_CONTAINER_INTERFACE_TOML_HH
#define GENERIC_CONTAINER_INTERFACE_TOML_HH

#ifdef __clang__
#pragma clang diagnostic ignored "-Wpoison-system-directories"
#endif

#include "GenericContainer.hh"

#include <fstream>

namespace GC_namespace
{

  using std::ifstream;

  //!
  //! \addtogroup TOML
  //!
  //! @{

  //!
  //! Convert a TOML file  to a `GenericContainer`.
  //!
  //! This function reads TOML data from the provided input stream and
  //! populates the given `GenericContainer` with the parsed data.
  //!
  //! \param[in]  file_name Input file name containing TOML data.
  //! \param[out] gc        The `GenericContainer` to be populated.
  //! \return     true      if the conversion was successful, false otherwise.
  //!
  inline bool file_TOML_to_GC( string_view file_name, GenericContainer & gc )
  {
    ifstream stream( file_name.data() );
    gc.clear();
    return gc.from_toml( stream );
  }

  //!
  //! Convert a TOML file stream to a `GenericContainer`.
  //!
  //! This function reads TOML data from the provided input stream and
  //! populates the given `GenericContainer` with the parsed data.
  //!
  //! \param[in]  stream Input stream containing TOML data.
  //! \param[out] gc     The `GenericContainer` to be populated.
  //! \return     true   if the conversion was successful, false otherwise.
  //!
  inline bool TOML_to_GC( istream_type & stream, GenericContainer & gc )
  {
    gc.clear();
    return gc.from_toml( stream );
  }

  //!
  //! Convert a TOML string to a `GenericContainer`.
  //!
  //! This function parses the given TOML string and populates the
  //! specified `GenericContainer` with the resulting data.
  //!
  //! \param[in]  DATA The TOML string to parse.
  //! \param[out] gc   The `GenericContainer` to be populated.
  //! \return     true if the conversion was successful, false otherwise.
  //!
  inline bool TOML_to_GC( string const & DATA, GenericContainer & gc )
  {
    istringstream stream( DATA );
    gc.clear();
    return gc.from_toml( stream );
  }

  //!
  //! Convert a vector of TOML strings to a `GenericContainer`.
  //!
  //! This function parses each string in the given vector as TOML
  //! and populates the specified `GenericContainer` with the resulting
  //! data.
  //!
  //! \param[in]  gc   The `GenericContainer` to convert.
  //! \param[out] DATA String to store the TOML encoded GenericContainer.
  //!
  inline void GC_to_TOML( GenericContainer const & gc, std::string & res )
  {
    ostringstream stream;
    gc.to_toml( stream );
    res = stream.str();
  }

  //!
  //! Convert a `GenericContainer` to a TOML file stream.
  //!
  //! This function converts the contents of the provided `GenericContainer`
  //! into TOML format and writes it to the specified output stream.
  //!
  //! \param[in]  gc     The `GenericContainer` to convert.
  //! \param[out] stream Output stream to write the TOML data.
  //!
  inline void GC_to_TOML( GenericContainer const & gc, ostream_type & stream )
  {
    gc.to_toml( stream );
  }

  //!
  //! @}
  //!
}  // namespace GC_namespace

#endif
