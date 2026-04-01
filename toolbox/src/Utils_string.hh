/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2026                                                      |
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
// file: Utils_string.hh
//
// \brief UTF-8 string manipulation utilities and border drawing characters
// \author Enrico Bertolazzi
// \date 2017
// \note Revised 2026 - Simplified and improved UTF-8 width calculation
//

#pragma once
#ifndef UTILS_STRING_HH
#define UTILS_STRING_HH

#include "Utils.hh"

namespace Utils
{

  // ============================================================================
  // STRING CASE UTILITIES
  // ============================================================================

  /**
   * \brief Converts a string to uppercase in-place
   * \param s String to convert (modified in place)
   * \note Uses std::toupper for each character
   */
  inline void to_upper( string & s )
  {
    std::transform( s.begin(), s.end(), s.begin(), []( unsigned char c ) { return std::toupper( c ); } );
  }

  /**
   * \brief Converts a string to lowercase in-place
   * \param s String to convert (modified in place)
   * \note Uses std::tolower for each character
   */
  inline void to_lower( string & s )
  {
    std::transform( s.begin(), s.end(), s.begin(), []( unsigned char c ) { return std::tolower( c ); } );
  }

  /**
   * \brief Checks if all characters in a string are lowercase
   * \param s String to check
   * \return true if all characters are lowercase, false otherwise
   */
  inline bool is_lower( string_view s )
  {
    return std::all_of( s.begin(), s.end(), islower );
  }

  /**
   * \brief Checks if all characters in a string are uppercase
   * \param s String to check
   * \return true if all characters are uppercase, false otherwise
   */
  inline bool is_upper( string_view s )
  {
    return std::all_of( s.begin(), s.end(), isupper );
  }

  /**
   * \brief Checks if all characters in a string are alphabetic
   * \param s String to check
   * \return true if all characters are alphabetic, false otherwise
   */
  inline bool is_alpha( string_view s )
  {
    return std::all_of( s.begin(), s.end(), isalpha );
  }

  /**
   * \brief Checks if all characters in a string are alphanumeric
   * \param s String to check
   * \return true if all characters are alphanumeric, false otherwise
   */
  inline bool is_alphanum( string_view s )
  {
    return std::all_of( s.begin(), s.end(), isalnum );
  }

  /**
   * \brief Checks if all characters in a string are decimal digits
   * \param s String to check
   * \return true if all characters are digits, false otherwise
   */
  inline bool is_digits( string_view s )
  {
    return std::all_of( s.begin(), s.end(), isdigit );
  }

  /**
   * \brief Checks if all characters in a string are hexadecimal digits
   * \param s String to check
   * \return true if all characters are hexadecimal digits, false otherwise
   */
  inline bool is_xdigits( string_view s )
  {
    return std::all_of( s.begin(), s.end(), isxdigit );
  }

  // ============================================================================
  // UTF-8 LOW LEVEL HELPERS
  // ============================================================================

  /**
   * \brief Determines the length in bytes of a UTF-8 character from its first byte
   * \param c First byte of UTF-8 character
   * \return Number of bytes in UTF-8 sequence (1-4), 0 for invalid start byte
   */
  inline int utf8_char_length( unsigned char c )
  {
    if ( c < 0x80 ) return 1;
    if ( ( c & 0xE0 ) == 0xC0 ) return 2;
    if ( ( c & 0xF0 ) == 0xE0 ) return 3;
    if ( ( c & 0xF8 ) == 0xF0 ) return 4;
    return 0;
  }

  /**
   * \brief Checks if a byte is a UTF-8 continuation byte
   * \param c Byte to check
   * \return true if byte is a continuation byte (binary pattern 10xxxxxx)
   */
  inline bool utf8_is_continuation( unsigned char c )
  {
    return ( c & 0xC0 ) == 0x80;
  }

  /**
   * \brief Decodes the next UTF-8 codepoint from a string
   * \param s UTF-8 encoded string
   * \param i Index in string (input/output, will be advanced past the character)
   * \return Unicode codepoint, or -1 on error (invalid UTF-8 sequence)
   */
  inline int utf8_next( std::string const & s, int & i )
  {
    int ss = static_cast<int>( s.size() );
    if ( i >= ss ) return -1;

    unsigned char c = static_cast<unsigned char>( s[i] );

    if ( c < 0x80 )
    {
      ++i;
      return c;
    }

    int len = utf8_char_length( c );
    if ( len == 0 || i + len > ss )
    {
      ++i;
      return -1;
    }

    for ( int k = 1; k < len; ++k )
      if ( !utf8_is_continuation( static_cast<unsigned char>( s[i + k] ) ) )
      {
        ++i;
        return -1;
      }

    int cp = 0;
    if ( len == 2 )
      cp = ( ( c & 0x1F ) << 6 ) | ( s[i + 1] & 0x3F );
    else if ( len == 3 )
      cp = ( ( c & 0x0F ) << 12 ) | ( ( s[i + 1] & 0x3F ) << 6 ) | ( s[i + 2] & 0x3F );
    else
      cp = ( ( c & 0x07 ) << 18 ) | ( ( s[i + 1] & 0x3F ) << 12 ) | ( ( s[i + 2] & 0x3F ) << 6 ) | ( s[i + 3] & 0x3F );

    i += len;
    return cp;
  }

  // ============================================================================
  // UNICODE WIDTH
  // ============================================================================

  /**
   * \brief Calculates display width (terminal columns) of a Unicode codepoint
   * \param cp Unicode codepoint
   * \return 0, 1 or 2 columns based on character category and rendering rules
   * \note This is a simplified implementation; full East Asian Width classification
   *       would require more complex tables.
   */
  inline int utf8_character_width( int cp )
  {
    if ( cp < 0 ) return 0;

    // Control characters (C0, C1, DELETE)
    if ( cp < 0x20 || ( cp >= 0x7F && cp < 0xA0 ) ) return 0;

    // Combining marks (zero width)
    if (
      ( cp >= 0x0300 && cp <= 0x036F ) || ( cp >= 0x1AB0 && cp <= 0x1AFF ) || ( cp >= 0x1DC0 && cp <= 0x1DFF ) ||
      ( cp >= 0x20D0 && cp <= 0x20FF ) || ( cp >= 0xFE20 && cp <= 0xFE2F ) )
      return 0;

    // ---- Explicit emoji/symbols shown as wide in modern terminals ----
    if (
      cp == 0x2705 ||                      // ✅ White Heavy Check Mark
      cp == 0x274C ||                      // ❌ Cross Mark - ADDED
      cp == 0x26A1 ||                      // ⚡ High Voltage Sign
      ( cp >= 0x2795 && cp <= 0x2797 ) ||  // ➕➖➗ (heavy plus/minus/divide)
      ( cp >= 0x2B05 && cp <= 0x2B07 ) ||  // ⬅⬇⬆ (heavy arrows)
      cp == 0x2B50 ||                      // ⭐ White Medium Star
      cp == 0x2B55                         // ⭕ Heavy Large Circle
    )
      return 2;

    // CJK / Fullwidth characters (East Asian Width property "W" or "F")
    if (
      ( cp >= 0x1100 && cp <= 0x115F ) || ( cp >= 0x2329 && cp <= 0x232A ) ||
      ( cp >= 0x2E80 && cp <= 0xA4CF && cp != 0x303F ) || ( cp >= 0xAC00 && cp <= 0xD7AF ) ||
      ( cp >= 0xF900 && cp <= 0xFAFF ) || ( cp >= 0xFE10 && cp <= 0xFE6F ) || ( cp >= 0xFF00 && cp <= 0xFF60 ) ||
      ( cp >= 0xFFE0 && cp <= 0xFFE6 ) )
      return 2;

    // Emoji blocks (Miscellaneous Symbols and Pictographs, etc.)
    if ( cp >= 0x1F300 && cp <= 0x1F9FF ) return 2;

    // Supplementary planes (typically wide characters)
    if ( cp >= 0x20000 ) return 2;

    return 1;
  }

  /**
   * \brief Repeats a string n times
   * \param s String to repeat
   * \param n Number of repetitions
   * \return Concatenated string with n repetitions of s
   */
  inline std::string repeat( const std::string & s, int n )
  {
    std::string r;
    r.reserve( s.length() * n );
    for ( int i = 0; i < n; ++i ) r += s;
    return r;
  }

  // ============================================================================
  // UTF-8 STRING UTILITIES
  // ============================================================================

  /**
   * \brief Calculates display width of a UTF-8 string in terminal columns
   * \param s UTF-8 encoded string
   * \return Total display width (sum of widths of individual Unicode characters)
   */
  inline int utf8_display_width( std::string const & s )
  {
    int w  = 0;
    int ss = static_cast<int>( s.size() );
    for ( int i = 0; i < ss; ) w += utf8_character_width( utf8_next( s, i ) );
    return w;
  }

  /**
   * \brief Truncates a UTF-8 string to fit within specified display width
   * \param s UTF-8 encoded string
   * \param max_width Maximum allowed display width
   * \param ellipsis String appended when truncated (default empty)
   * \return Truncated string (with ellipsis if needed)
   * \note The ellipsis width is counted in the total; the result width ≤ max_width.
   */
  inline std::string utf8_truncate( std::string const & s, int max_width, std::string const & ellipsis = "" )
  {
    std::string out;
    int         i = 0;
    int         w = 0;

    int ew = utf8_display_width( ellipsis );
    int ss = static_cast<int>( s.size() );

    while ( i < ss )
    {
      int prev = i;
      int cp   = utf8_next( s, i );
      int cw   = utf8_character_width( cp );

      if ( w + cw + ew > max_width ) break;

      out.append( s, prev, i - prev );
      w += cw;
    }

    if ( i < ss ) out += ellipsis;
    return out;
  }

  /**
   * \brief Pads a UTF-8 string to reach specified display width
   * \param s UTF-8 encoded string
   * \param width Desired total display width
   * \param pad Padding string (default single space)
   * \return Padded string (original + padding up to width)
   * \note If current width ≥ width, returns original string.
   */
  inline std::string utf8_padding( std::string const & s, int width, std::string const & pad = " " )
  {
    int cur = utf8_display_width( s );
    if ( cur >= width || pad.empty() ) return s;

    int pw = utf8_display_width( pad );
    if ( pw <= 0 ) return s;

    std::string out = s;
    while ( cur + pw <= width )
    {
      out += pad;
      cur += pw;
    }
    return out;
  }

  // ============================================================================
  // BOX DRAWING CHARACTERS - SINGLE LINE
  // ============================================================================

  constexpr char const * border_top       = "─"; /**< Top border (horizontal) */
  constexpr char const * border_top_mid   = "┬"; /**< Top-middle junction (T-down) */
  constexpr char const * border_top_left  = "┌"; /**< Top-left corner */
  constexpr char const * border_top_right = "┐"; /**< Top-right corner */

  constexpr char const * border_bottom       = "─"; /**< Bottom border (horizontal) */
  constexpr char const * border_bottom_mid   = "┴"; /**< Bottom-middle junction (T-up) */
  constexpr char const * border_bottom_left  = "└"; /**< Bottom-left corner */
  constexpr char const * border_bottom_right = "┘"; /**< Bottom-right corner */

  constexpr char const * border_left     = "│"; /**< Left border (vertical) */
  constexpr char const * border_left_mid = "├"; /**< Left-middle junction (T-right) */

  constexpr char const * border_mid     = "─"; /**< Middle border (horizontal) */
  constexpr char const * border_mid_mid = "┼"; /**< Center junction (cross) */

  constexpr char const * border_right     = "│"; /**< Right border (vertical) */
  constexpr char const * border_right_mid = "┤"; /**< Right-middle junction (T-left) */

  constexpr char const * border_middle = "│"; /**< Vertical separator between columns */

  // ============================================================================
  // BOX DRAWING CHARACTERS - BOLD/THICK LINE
  // ============================================================================

  constexpr char const * border_top_bold       = "━"; /**< Bold top border */
  constexpr char const * border_top_mid_bold   = "┳"; /**< Bold top-middle junction */
  constexpr char const * border_top_left_bold  = "┏"; /**< Bold top-left corner */
  constexpr char const * border_top_right_bold = "┓"; /**< Bold top-right corner */

  constexpr char const * border_bottom_bold       = "━"; /**< Bold bottom border */
  constexpr char const * border_bottom_mid_bold   = "┻"; /**< Bold bottom-middle junction */
  constexpr char const * border_bottom_left_bold  = "┗"; /**< Bold bottom-left corner */
  constexpr char const * border_bottom_right_bold = "┛"; /**< Bold bottom-right corner */

  constexpr char const * border_left_bold     = "┃"; /**< Bold left border */
  constexpr char const * border_left_mid_bold = "┣"; /**< Bold left-middle junction */

  constexpr char const * border_mid_bold     = "━"; /**< Bold middle border */
  constexpr char const * border_mid_mid_bold = "╋"; /**< Bold center junction */

  constexpr char const * border_right_bold     = "┃"; /**< Bold right border */
  constexpr char const * border_right_mid_bold = "┫"; /**< Bold right-middle junction */

  constexpr char const * border_middle_bold = "┃"; /**< Bold vertical separator */

  // ============================================================================
  // BOX DRAWING CHARACTERS - DOUBLE LINE
  // ============================================================================

  constexpr char const * border_top_double       = "═"; /**< Double top border */
  constexpr char const * border_top_mid_double   = "╦"; /**< Double top-middle junction */
  constexpr char const * border_top_left_double  = "╔"; /**< Double top-left corner */
  constexpr char const * border_top_right_double = "╗"; /**< Double top-right corner */

  constexpr char const * border_bottom_double       = "═"; /**< Double bottom border */
  constexpr char const * border_bottom_mid_double   = "╩"; /**< Double bottom-middle junction */
  constexpr char const * border_bottom_left_double  = "╚"; /**< Double bottom-left corner */
  constexpr char const * border_bottom_right_double = "╝"; /**< Double bottom-right corner */

  constexpr char const * border_left_double     = "║"; /**< Double left border */
  constexpr char const * border_left_mid_double = "╠"; /**< Double left-middle junction */

  constexpr char const * border_mid_double     = "═"; /**< Double middle border */
  constexpr char const * border_mid_mid_double = "╬"; /**< Double center junction */

  constexpr char const * border_right_double     = "║"; /**< Double right border */
  constexpr char const * border_right_mid_double = "╣"; /**< Double right-middle junction */

  constexpr char const * border_middle_double = "║"; /**< Double vertical separator */


}  // namespace Utils

#endif

//
// eof: Utils_string.hh
//
