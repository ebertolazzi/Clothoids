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

#ifndef GENERIC_CONTAINER_INTERFACE_YAML_HH
#define GENERIC_CONTAINER_INTERFACE_YAML_HH

#ifdef __clang__
#pragma clang diagnostic ignored "-Wpoison-system-directories"
#endif

#include "GenericContainer.hh"

#include <fstream>

namespace GC_namespace {

  using std::ifstream;

  //!
  //! \addtogroup YAML
  //!
  //! @{

  //!
  //! Convert a YAML file  to a `GenericContainer`.
  //!
  //! This function reads YAML data from the provided input stream and
  //! populates the given `GenericContainer` with the parsed data.
  //!
  //! \param[in]  file_name Input file name containing YAML data.
  //! \param[out] gc        The `GenericContainer` to be populated.
  //! \return     true      if the conversion was successful, false otherwise.
  //!
  inline
  bool
  file_YAML_to_GC(
    string const     & file_name,
    GenericContainer & gc
  ) {
    ifstream stream(file_name.c_str());
    gc.clear();
    return gc.from_yaml( stream );
  }

  //!
  //! Convert a YAML file stream to a `GenericContainer`.
  //!
  //! This function reads YAML data from the provided input stream and
  //! populates the given `GenericContainer` with the parsed data.
  //!
  //! \param[in]  stream Input stream containing YAML data.
  //! \param[out] gc     The `GenericContainer` to be populated.
  //! \return     true   if the conversion was successful, false otherwise.
  //!
  inline
  bool
  YAML_to_GC(
    istream_type     & stream,
    GenericContainer & gc
  ) {
    gc.clear();
    return gc.from_yaml( stream );
  }

  //!
  //! Convert a YAML string to a `GenericContainer`.
  //!
  //! This function parses the given YAML string and populates the
  //! specified `GenericContainer` with the resulting data.
  //!
  //! \param[in]  DATA The YAML string to parse.
  //! \param[out] gc   The `GenericContainer` to be populated.
  //! \return     true if the conversion was successful, false otherwise.
  //!
  inline
  bool
  YAML_to_GC(
    string const     & DATA,
    GenericContainer & gc
  ) {
    istringstream stream(DATA);
    gc.clear();
    return gc.from_yaml(stream);
  }

  //!
  //! Convert a vector of YAML strings to a `GenericContainer`.
  //!
  //! This function parses each string in the given vector as YAML
  //! and populates the specified `GenericContainer` with the resulting
  //! data.
  //!
  //! \param[in]  gc   The `GenericContainer` to convert.
  //! \param[out] DATA String to store the YAML encoded GenericContainer.
  //!
  inline
  void
  GC_to_YAML(
    GenericContainer const & gc,
    string                 & DATA
  ) {
    ostringstream stream;
    gc.to_yaml( stream );
    DATA = stream.str();
  }

  //!
  //! Convert a `GenericContainer` to a YAML file stream.
  //!
  //! This function converts the contents of the provided `GenericContainer`
  //! into YAML format and writes it to the specified output stream.
  //!
  //! \param[in]  gc     The `GenericContainer` to convert.
  //! \param[out] stream Output stream to write the YAML data.
  //!
  inline
  void
  GC_to_YAML(
    GenericContainer const & gc,
    ostream_type           & stream
  ) {
    gc.to_yaml( stream );
  }

  //!
  //! @}
  //!
}

#endif
