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

namespace Utils
{

  inline string fmt_table_left_top()
  {
    return "┌";
  };
  inline string fmt_table_left_middle()
  {
    return "├";
  };
  inline string fmt_table_left_bottom()
  {
    return "└";
  };

  inline string fmt_table_middle_top()
  {
    return "┬";
  };
  inline string fmt_table_middle_middle()
  {
    return "┼";
  };
  inline string fmt_table_middle_bottom()
  {
    return "┴";
  };

  inline string fmt_table_right_top()
  {
    return "┐";
  };
  inline string fmt_table_right_middle()
  {
    return "┤";
  };
  inline string fmt_table_right_bottom()
  {
    return "┘";
  };

  inline string fmt_table_cross()
  {
    return "┼";
  };

  inline string fmt_table_bar()
  {
    return "─";
  }
  inline string fmt_table_vbar()
  {
    return "│";
  }
  inline string fmt_table_dot()
  {
    return "•";
  }
  inline string fmt_table_vdots()
  {
    return "⋮";
  }

  using std::string_view;

  string fmt_table_row(
    unsigned    width,
    string_view L,
    string_view R,
    string_view F,
    string_view title,
    string_view align );
  string fmt_table_row(
    unsigned                           width,
    string_view                        L,
    string_view                        C,
    string_view                        R,
    string_view                        F,
    std::initializer_list<string_view> names,
    string_view                        align );
  string fmt_table_row( unsigned width, string_view L, string_view C, string_view R, string_view F, unsigned N );

  inline string fmt_table_row( unsigned width, string_view title = "", string_view align = "^", string_view fill = " " )
  {
    return fmt_table_row( width, "│", "│\n", fill, title, align );
  }
  inline string fmt_table_top_row(
    unsigned    width,
    string_view title = "",
    string_view align = "^",
    string_view fill  = "─" )
  {
    return fmt_table_row( width, "┌", "┐\n", fill, title, align );
  }
  inline string fmt_table_middle_row(
    unsigned    width,
    string_view title = "",
    string_view align = "^",
    string_view fill  = "─" )
  {
    return fmt_table_row( width, "├", "┤\n", fill, title, align );
  }
  inline string fmt_table_bottom_row(
    unsigned    width,
    string_view title = "",
    string_view align = "^",
    string_view fill  = "─" )
  {
    return fmt_table_row( width, "└", "┘\n", fill, title, align );
  }

  inline string fmt_table_row(
    unsigned                           width,
    std::initializer_list<string_view> names,
    string_view                        align = "<",
    string_view                        fill  = " " )
  {
    return fmt_table_row( width, "│", "│", "│\n", fill, names, align );
  }
  inline string fmt_table_top_row(
    unsigned                           width,
    std::initializer_list<string_view> names,
    string_view                        align = "<",
    string_view                        fill  = "─" )
  {
    return fmt_table_row( width, "┌", "─", "┐\n", fill, names, align );
  }
  inline string fmt_table_middle_row(
    unsigned                           width,
    std::initializer_list<string_view> names,
    string_view                        align = "<",
    string_view                        fill  = "─" )
  {
    return fmt_table_row( width, "├", "┼", "┤\n", fill, names, align );
  }
  inline string fmt_table_bottom_row(
    unsigned                           width,
    std::initializer_list<string_view> names,
    string_view                        align = "<",
    string_view                        fill  = "─" )
  {
    return fmt_table_row( width, "└", "─", "┘\n", fill, names, align );
  }

  inline string fmt_table_row( unsigned width, unsigned N, string_view fill = " " )
  {
    return fmt_table_row( width, "│", "│", "│\n", fill, N );
  }
  inline string fmt_table_top_row( unsigned width, unsigned N, string_view fill = "─" )
  {
    return fmt_table_row( width, "┌", "┬", "┐\n", fill, N );
  }
  inline string fmt_table_middle_row( unsigned width, unsigned N, string_view fill = "─" )
  {
    return fmt_table_row( width, "├", "┼", "┤\n", fill, N );
  }
  inline string fmt_table_bottom_row( unsigned width, unsigned N, string_view fill = "─" )
  {
    return fmt_table_row( width, "└", "┴", "┘\n", fill, N );
  }

  // Helper for vector formatting
  template <typename Scalar> inline string format_index_vector( std::vector<Scalar> const & v, size_t max_size = 20 )
  {
    std::string tmp{ "[" };
    size_t      v_size = v.size();
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

  /**
   * @brief Format a vector of indices in a compact representation
   *
   * @tparam T Index type (typically size_t or int)
   * @param indices Vector of indices to format
   * @param max_display Maximum number of indices to display before truncating
   * @return std::string Compact string representation
   */
  template <typename T>
  inline std::string format_index_vector_compact( std::vector<T> const & indices, size_t max_display = 5 )
  {
    if ( indices.empty() ) return "[]";

    std::stringstream ss;
    ss << "[";

    if ( indices.size() <= max_display )
    {
      for ( size_t i = 0; i < indices.size(); ++i )
      {
        if ( i > 0 ) ss << ", ";
        ss << indices[i];
      }
    }
    else
    {
      for ( size_t i = 0; i < max_display - 1; ++i )
      {
        if ( i > 0 ) ss << ", ";
        ss << indices[i];
      }
      ss << ", ..., " << indices.back();
    }

    ss << "]";
    return ss.str();
  }

  template <typename Vector> inline std::string format_reduced_vector( Vector const & v, size_t max_size = 10 )
  {
    std::string tmp{ "[" };
    size_t      v_size = v.size();

    if ( v_size == 0 ) { return "[]"; }

    if ( v_size <= max_size )
    {
      for ( size_t i = 0; i < v_size; ++i ) tmp += fmt::format( "{:.4f}, ", v( i ) );
    }
    else
    {
      for ( size_t i = 0; i < max_size - 3; ++i ) tmp += fmt::format( "{:.4f}, ", v( i ) );
      tmp += "..., ";
      for ( size_t i = v_size - 3; i < v_size; ++i ) tmp += fmt::format( "{:.4f}, ", v( i ) );
    }

    if ( v_size > 0 )
    {
      tmp.pop_back();
      tmp.pop_back();
    }
    tmp += "]";
    return tmp;
  }

  namespace PrintColors
  {
    constexpr auto HEADER    = fmt::fg( fmt::color::light_blue );
    constexpr auto SUCCESS   = fmt::fg( fmt::color::green );
    constexpr auto WARNING   = fmt::fg( fmt::color::yellow );
    constexpr auto ERROR     = fmt::fg( fmt::color::red );
    constexpr auto INFO      = fmt::fg( fmt::color::cyan );
    constexpr auto ITERATION = fmt::fg( fmt::color::white );
    constexpr auto DETAIL    = fmt::fg( fmt::color::gray );
  }  // namespace PrintColors

}  // namespace Utils

#endif

//
// EOF: Utils_fmt.hh
//
