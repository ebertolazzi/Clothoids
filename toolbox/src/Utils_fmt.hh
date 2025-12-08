/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2025                                                      |
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
// file: Utils_fmt.hh
//

#pragma once

#ifndef UTILS_FMT_dot_HH
#define UTILS_FMT_dot_HH

#if defined( __llvm__ ) || defined( __clang__ )
#pragma clang diagnostic ignored "-Wduplicate-enum"
#endif

#include "Utils.hh"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "Utils/3rd/spdlog/spdlog.h"
#include "Utils/3rd/spdlog/fmt/bundled/std.h"
#include "Utils/3rd/spdlog/fmt/bundled/chrono.h"
#include "Utils/3rd/spdlog/fmt/bundled/color.h"
#include "Utils/3rd/spdlog/fmt/bundled/ostream.h"
#include "Utils/3rd/spdlog/fmt/bundled/printf.h"
#endif

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include <string.h>
#ifndef __FILENAME__
#define __FILENAME__ ( strrchr( __FILE__, '/' ) ? strrchr( "/" __FILE__, '/' ) + 1 : __FILE__ )
#endif

#ifndef UTILS_ERROR0
#define UTILS_ERROR0( MSG ) throw Utils::Runtime_Error( MSG, __FILENAME__, __LINE__ )
#endif

#ifndef UTILS_ASSERT0
#define UTILS_ASSERT0( COND, MSG )                                                                                     \
  if ( !( COND ) ) UTILS_ERROR0( MSG )
#endif

#ifndef UTILS_WARNING0
#define UTILS_WARNING0( COND, MSG )                                                                                    \
  if ( !( COND ) ) std::cerr << MSG
#endif

#ifndef UTILS_ERROR
#define UTILS_ERROR( ... ) throw Utils::Runtime_Error( fmt::format( __VA_ARGS__ ), __FILENAME__, __LINE__ )
#endif

#ifndef UTILS_ASSERT
#define UTILS_ASSERT( COND, ... )                                                                                      \
  if ( !( COND ) ) UTILS_ERROR( __VA_ARGS__ )
#endif

#ifndef UTILS_WARNING
#define UTILS_WARNING( COND, ... )                                                                                     \
  if ( !( COND ) ) fmt::print( __VA_ARGS__ )
#endif

#ifdef UTILS_NO_DEBUG
#ifndef UTILS_ASSERT0_DEBUG
#define UTILS_ASSERT0_DEBUG( COND, MSG )
#endif
#ifndef UTILS_ASSERT_DEBUG
#define UTILS_ASSERT_DEBUG( COND, ... )
#endif
#else
#ifndef UTILS_ASSERT0_DEBUG
#define UTILS_ASSERT0_DEBUG( COND, MSG ) UTILS_ASSERT0( COND, MSG )
#endif
#ifndef UTILS_ASSERT_DEBUG
#define UTILS_ASSERT_DEBUG( COND, ... ) UTILS_ASSERT( COND, __VA_ARGS__ )
#endif
#endif

#endif

#include <string>
#include <string_view>

namespace Utils
{

  using std::runtime_error;

  //!
  //! \brief Custom runtime error class for handling runtime exceptions.
  //!
  //! This class extends the standard `std::runtime_error` to include additional
  //! context information, specifically the file name and line number where the
  //! error occurred. It provides constructors that accept a reason for the
  //! error and formats the error message accordingly.
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

  class Runtime_Error : public runtime_error
  {
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
    explicit Runtime_Error( string_view reason, string_view file, int line )
      : std::runtime_error( fmt::format( "\n{}\nOn File:{}:{}\n", reason, file, line ) )
    {
    }

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

  inline string
  fmt_table_left_top()
  {
    return "┌";
  };
  inline string
  fmt_table_left_middle()
  {
    return "├";
  };
  inline string
  fmt_table_left_bottom()
  {
    return "└";
  };

  inline string
  fmt_table_middle_top()
  {
    return "┬";
  };
  inline string
  fmt_table_middle_middle()
  {
    return "┼";
  };
  inline string
  fmt_table_middle_bottom()
  {
    return "┴";
  };

  inline string
  fmt_table_right_top()
  {
    return "┐";
  };
  inline string
  fmt_table_right_middle()
  {
    return "┤";
  };
  inline string
  fmt_table_right_bottom()
  {
    return "┘";
  };

  inline string
  fmt_table_cross()
  {
    return "┼";
  };

  inline string
  fmt_table_bar()
  {
    return "─";
  }
  inline string
  fmt_table_vbar()
  {
    return "│";
  }
  inline string
  fmt_table_dot()
  {
    return "•";
  }
  inline string
  fmt_table_vdots()
  {
    return "⋮";
  }

  using std::string_view;

  string fmt_table_row( unsigned    width,
                        string_view L,
                        string_view R,
                        string_view F,
                        string_view title,
                        string_view align );
  string fmt_table_row( unsigned                           width,
                        string_view                        L,
                        string_view                        C,
                        string_view                        R,
                        string_view                        F,
                        std::initializer_list<string_view> names,
                        string_view                        align );
  string fmt_table_row( unsigned width, string_view L, string_view C, string_view R, string_view F, unsigned N );

  inline string
  fmt_table_row( unsigned width, string_view title = "", string_view align = "^", string_view fill = " " )
  {
    return fmt_table_row( width, "│", "│\n", fill, title, align );
  }
  inline string
  fmt_table_top_row( unsigned width, string_view title = "", string_view align = "^", string_view fill = "─" )
  {
    return fmt_table_row( width, "┌", "┐\n", fill, title, align );
  }
  inline string
  fmt_table_middle_row( unsigned width, string_view title = "", string_view align = "^", string_view fill = "─" )
  {
    return fmt_table_row( width, "├", "┤\n", fill, title, align );
  }
  inline string
  fmt_table_bottom_row( unsigned width, string_view title = "", string_view align = "^", string_view fill = "─" )
  {
    return fmt_table_row( width, "└", "┘\n", fill, title, align );
  }

  inline string
  fmt_table_row( unsigned                           width,
                 std::initializer_list<string_view> names,
                 string_view                        align = "<",
                 string_view                        fill  = " " )
  {
    return fmt_table_row( width, "│", "│", "│\n", fill, names, align );
  }
  inline string
  fmt_table_top_row( unsigned                           width,
                     std::initializer_list<string_view> names,
                     string_view                        align = "<",
                     string_view                        fill  = "─" )
  {
    return fmt_table_row( width, "┌", "─", "┐\n", fill, names, align );
  }
  inline string
  fmt_table_middle_row( unsigned                           width,
                        std::initializer_list<string_view> names,
                        string_view                        align = "<",
                        string_view                        fill  = "─" )
  {
    return fmt_table_row( width, "├", "┼", "┤\n", fill, names, align );
  }
  inline string
  fmt_table_bottom_row( unsigned                           width,
                        std::initializer_list<string_view> names,
                        string_view                        align = "<",
                        string_view                        fill  = "─" )
  {
    return fmt_table_row( width, "└", "─", "┘\n", fill, names, align );
  }

  inline string
  fmt_table_row( unsigned width, unsigned N, string_view fill = " " )
  {
    return fmt_table_row( width, "│", "│", "│\n", fill, N );
  }
  inline string
  fmt_table_top_row( unsigned width, unsigned N, string_view fill = "─" )
  {
    return fmt_table_row( width, "┌", "┬", "┐\n", fill, N );
  }
  inline string
  fmt_table_middle_row( unsigned width, unsigned N, string_view fill = "─" )
  {
    return fmt_table_row( width, "├", "┼", "┤\n", fill, N );
  }
  inline string
  fmt_table_bottom_row( unsigned width, unsigned N, string_view fill = "─" )
  {
    return fmt_table_row( width, "└", "┴", "┘\n", fill, N );
  }

  // Helper for vector formatting
  template <typename Scalar>
  inline string
  format_index_vector( vector<Scalar> const & v, size_t max_size = 20 )
  {
    string tmp{ "[" };
    size_t v_size = v.size();
    if ( v_size <= max_size )
    {
      for ( size_t i = 0; i < v_size; ++i ) tmp += fmt::format( "{}, ", v[i] );
    }
    else
    {
      for ( size_t i{ 0 }; i < max_size - 3; ++i ) tmp += fmt::format( "{}, ", v[i] );
      tmp.pop_back();
      tmp += "..., ";
      for ( size_t i{ v_size - 3 }; i < v_size; ++i ) tmp += fmt::format( "{}, ", v[i] );
    }
    tmp.pop_back();
    tmp.pop_back();
    tmp += "]";
    return tmp;
  }
}  // namespace Utils

#endif

//
// EOF: Utils_fmt.hh
//
