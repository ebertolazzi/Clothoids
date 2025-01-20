/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2022                                                      |
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
// file: Utils_fmt.hh
//

#pragma once

#ifndef UTILS_FMT_dot_HH
#define UTILS_FMT_dot_HH

#include "Utils.hh"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "Utils/fmt/printf.h"
#include "Utils/fmt/chrono.h"
#include "Utils/fmt/ostream.h"
#include "Utils/fmt/color.h"
#include "Utils/fmt/std.h"
#endif

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include <string.h>
#ifndef __FILENAME__
  #define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr("/" __FILE__, '/') + 1 : __FILE__)
#endif

#ifndef UTILS_ERROR0
  #define UTILS_ERROR0(MSG) \
  throw Utils::Runtime_Error( MSG, __FILENAME__, __LINE__ )
#endif

#ifndef UTILS_ASSERT0
  #define UTILS_ASSERT0(COND,MSG) if ( !(COND) ) UTILS_ERROR0( MSG )
#endif

#ifndef UTILS_WARNING0
  #define UTILS_WARNING0(COND,MSG) if ( !(COND) ) std::cerr << MSG
#endif

#ifndef UTILS_ERROR
  #define UTILS_ERROR(...) \
  throw Utils::Runtime_Error( fmt::format(__VA_ARGS__), __FILENAME__, __LINE__ )
#endif

#ifndef UTILS_ASSERT
  #define UTILS_ASSERT(COND,...) if ( !(COND) ) UTILS_ERROR( __VA_ARGS__ )
#endif

#ifndef UTILS_WARNING
  #define UTILS_WARNING(COND,...) if ( !(COND) ) fmt::print( __VA_ARGS__ )
#endif

#ifdef UTILS_NO_DEBUG
  #ifndef UTILS_ASSERT0_DEBUG
    #define UTILS_ASSERT0_DEBUG(COND,MSG)
  #endif
  #ifndef UTILS_ASSERT_DEBUG
    #define UTILS_ASSERT_DEBUG(COND,...)
  #endif
#else
  #ifndef UTILS_ASSERT0_DEBUG
    #define UTILS_ASSERT0_DEBUG(COND,MSG) UTILS_ASSERT0(COND,MSG)
  #endif
  #ifndef UTILS_ASSERT_DEBUG
    #define UTILS_ASSERT_DEBUG(COND,...) UTILS_ASSERT(COND,__VA_ARGS__)
  #endif
#endif

#endif

namespace Utils {

  using std::runtime_error;

  //!
  //! \brief Custom runtime error class for handling runtime exceptions.
  //!
  //! This class extends the standard `std::runtime_error` to include additional
  //! context information, specifically the file name and line number where the
  //! error occurred. It provides constructors that accept a reason for the error
  //! and formats the error message accordingly.
  //!
  //! **Usage**
  //!
  //! \code
  //! try {
  //!     throw Runtime_Error("An error occurred", __FILE__, __LINE__);
  //! } catch (const Runtime_Error& e) {
  //!     std::cerr << e.what();
  //! }
  //! \endcode

  class Runtime_Error : public runtime_error {
  public:
    //!
    //! \brief Constructs a Runtime_Error instance with a given reason.
    //!
    //! This constructor initializes the error with a specified reason,
    //! the file where the error occurred, and the line number. It formats
    //! the error message accordingly.
    //!
    //! \param reason A string that describes the reason for the error.
    //! \param file The name of the file where the error occurred.
    //! \param line The line number in the file where the error occurred.
    //!
    explicit
    Runtime_Error(
      std::string const & reason,
      char const *        file,
      int                 line
    )
    : std::runtime_error( fmt::format( "\n{}\nOn File:{}:{}\n", reason, file, line ) )
    { }

    //!
    //! \brief Constructs a Runtime_Error instance with a given reason.
    //!
    //! This constructor initializes the error with a specified reason,
    //! the file where the error occurred, and the line number. It formats
    //! the error message accordingly.
    //!
    //! \param reason A C-style string that describes the reason for the error.
    //! \param file The name of the file where the error occurred.
    //! \param line The line number in the file where the error occurred.
    //!
    explicit
    Runtime_Error(
      char const * reason,
      char const * file,
      int          line
    )
    : std::runtime_error( fmt::format( "\n{}\nOn File:{}:{}\n", reason, file, line ) )
    { }

    //!
    //! \brief Returns a C-style string describing the error.
    //!
    //! This method overrides the `what()` method from `std::runtime_error`
    //! to provide a more detailed error message, including the reason for
    //! the error, the file name, and the line number.
    //!
    //! \return A C-style string representing the error message.
    //!
    char const * what() const noexcept override;
  };

}

#endif

//
// EOF: Utils_fmt.hh
//
