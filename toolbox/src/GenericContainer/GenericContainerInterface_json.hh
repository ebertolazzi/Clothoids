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

#pragma once

#ifndef GENERIC_CONTAINER_INTERFACE_JSON_HH
#define GENERIC_CONTAINER_INTERFACE_JSON_HH

#include "GenericContainer.hh"

#include <fstream>

namespace GC_namespace {

  using std::ifstream;

  //!
  //! \addtogroup JSON
  //!
  //! @{

  //!
  //! Convert a JSON file to a `GenericContainer`.
  //!
  //! This function reads JSON data from the provided input stream and
  //! populates the given `GenericContainer` with the parsed data.
  //!
  //! \param[in]  file_name Input file name containing JSON data.
  //! \param[out] gc        The `GenericContainer` to be populated.
  //! \return     true      if the conversion was successful, false otherwise.
  //!
  inline
  bool
  file_JSON_to_GC(
    string const     & file_name,
    GenericContainer & gc
  ) {
    ifstream stream(file_name.c_str());
    gc.clear();
    return gc.from_json( stream );
  }

  //!
  //! Convert a JSON file stream to a `GenericContainer`.
  //!
  //! This function reads JSON data from the provided input stream and
  //! populates the given `GenericContainer` with the parsed data.
  //!
  //! \param[in]  stream Input stream containing JSON data.
  //! \param[out] gc     The `GenericContainer` to be populated.
  //! \return    true    if the conversion was successful, false otherwise.
  //!
  inline
  bool
  JSON_to_GC(
    istream_type     & stream,
    GenericContainer & gc
  ) {
    gc.clear();
    return gc.from_json( stream );
  }

  //!
  //! Convert a JSON string to a `GenericContainer`.
  //!
  //! This function parses the given JSON string and populates the
  //! specified `GenericContainer` with the resulting data.
  //!
  //! \param[in]  DATA The JSON string to parse.
  //! \param[out] gc   The `GenericContainer` to be populated.
  //! \return     true if the conversion was successful, false otherwise.
  //!
  inline
  bool
  JSON_to_GC(
    string const     & DATA,
    GenericContainer & gc
  ) {
    istringstream stream(DATA);
    gc.clear();
    return gc.from_json(stream);
  }

  //!
  //! Convert a `GenericContainer` to a JSON string.
  //!
  //! This function converts the contents of the provided `GenericContainer`
  //! into JSON format and stores it in the specified string.
  //!
  //! \param[in]  gc   The `GenericContainer` to convert.
  //! \param[out] DATA String to store the JSON encoded GenericContainer.
  //!
  inline
  void
  GC_to_JSON(
    GenericContainer const & gc,
    string                 & DATA
  ) {
    ostringstream stream(DATA);
    gc.to_json( stream );
  }

  //!
  //! Convert a `GenericContainer` to a JSON file stream.
  //!
  //! This function converts the contents of the provided `GenericContainer`
  //! into JSON format and writes it to the specified output stream.
  //!
  //! \param[in]  gc     The `GenericContainer` to convert.
  //! \param[out] stream Output stream to write the JSON data.
  //!
  inline
  void
  GC_to_JSON(
    GenericContainer const & gc,
    ostream_type           & stream
  ) {
    gc.to_json( stream );
  }

  //!
  //! @}
  //!

}

#endif

