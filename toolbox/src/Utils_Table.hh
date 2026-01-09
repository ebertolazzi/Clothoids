/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2021                                                      |
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
 |      Via Sommarive 9, I-38123 Povo, Trento, Italy                        |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

//
// file: Utils_Table.hh
//

#pragma once
#ifndef UTILS_TABLE_HH
#define UTILS_TABLE_HH

#include "Utils.hh"
#include "Utils_string.hh"

/*\

  Based on terminal-table:

  https://github.com/Bornageek/terminal-table

  Copyright 2015 Andreas Wilhelm

  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

\*/

namespace Utils
{
  using std::string;
  using std::string_view;
  using std::vector;

  namespace Table
  {

    using integer = int;

    //!
    //! \brief Enum class defining alignment types for table cells.
    //!
    //! Provides options for aligning content inside table cells to the LEFT,
    //! RIGHT, or CENTER.
    //!
    using Alignment = enum class Table_align : integer { LEFT, RIGHT, CENTER };

    // =============================================================
    // Junction resolution
    // =============================================================
    struct Junction
    {
      bool up    = false;
      bool down  = false;
      bool left  = false;
      bool right = false;
    };

    inline std::string resolve_junction( Junction const & j )
    {
      if ( j.up && j.down && j.left && j.right ) return "┼";
      if ( j.up && j.down && j.right ) return "├";
      if ( j.up && j.down && j.left ) return "┤";
      if ( j.left && j.right && j.down ) return "┬";
      if ( j.left && j.right && j.up ) return "┴";
      if ( j.down && j.right ) return "┌";
      if ( j.down && j.left ) return "┐";
      if ( j.up && j.right ) return "└";
      if ( j.up && j.left ) return "┘";
      if ( j.left || j.right ) return "─";
      return "│";
    }

    //!
    //! \brief Defines the style and structure of table borders, padding, and
    //! alignment.
    //!
    //! The `Style` class configures visual properties of the table, including
    //! border characters (for all edges and dividers), cell padding, and
    //! alignment of text.
    //!
    class Style
    {
    private:
      string m_border_top       = "─";
      string m_border_top_mid   = "┬";
      string m_border_top_left  = "┌";
      string m_border_top_right = "┐";

      string m_border_bottom       = "─";
      string m_border_bottom_mid   = "┴";
      string m_border_bottom_left  = "└";
      string m_border_bottom_right = "┘";

      string m_border_left     = "│";
      string m_border_left_mid = "├";

      string m_border_mid     = "─";
      string m_border_mid_mid = "┼";

      string m_border_right     = "│";
      string m_border_right_mid = "┤";

      string m_border_middle = "│";

      integer m_padding_left  = 1;
      integer m_padding_right = 1;

      Alignment m_Align = Alignment::LEFT;

      integer m_Width = 0;

    public:
      //!
      //! \brief Default constructor initializing the table style with default
      //! borders.
      //!
      Style() = default;

      string border_top() const { return m_border_top; }
      void   border_top( string borderStyle ) { m_border_top = borderStyle; }

      string border_top_mid() const { return m_border_top_mid; }
      void   border_top_mid( string borderStyle ) { m_border_top_mid = borderStyle; }

      string border_top_left() const { return m_border_top_left; }
      void   border_top_left( string borderStyle ) { m_border_top_left = borderStyle; }

      string border_top_right() const { return m_border_top_right; }
      void   border_top_right( string borderStyle ) { m_border_top_right = borderStyle; }

      string border_bottom() const { return m_border_bottom; }
      void   border_bottom( string borderStyle ) { m_border_bottom = borderStyle; }

      string border_bottom_mid() const { return m_border_bottom_mid; }
      void   border_bottom_mid( string borderStyle ) { m_border_bottom_mid = borderStyle; }

      string border_bottom_left() const { return m_border_bottom_left; }
      void   border_bottom_left( string borderStyle ) { m_border_bottom_left = borderStyle; }

      string border_bottom_right() const { return m_border_bottom_right; }
      void   border_bottom_right( string borderStyle ) { m_border_bottom_right = borderStyle; }

      string border_left() const { return m_border_left; }
      void   border_left( string borderStyle ) { m_border_left = borderStyle; }

      string border_left_mid() const { return m_border_left_mid; }
      void   border_left_mid( string borderStyle ) { m_border_left_mid = borderStyle; }

      string border_mid() const { return m_border_mid; }
      void   border_mid( string borderStyle ) { m_border_mid = borderStyle; }

      string border_mid_mid() const { return m_border_mid_mid; }
      void   border_mid_mid( string borderStyle ) { m_border_mid_mid = borderStyle; }

      string border_right() const { return m_border_right; }
      void   border_right( string borderStyle ) { m_border_right = borderStyle; }

      string border_right_mid() const { return m_border_right_mid; }
      void   border_right_mid( string borderStyle ) { m_border_right_mid = borderStyle; }

      string border_middle() const { return m_border_middle; }
      void   border_middle( string borderStyle ) { m_border_middle = borderStyle; }

      integer padding_left() const { return m_padding_left; }
      void    padding_left( integer padding ) { m_padding_left = padding; }

      integer padding_right() const { return m_padding_right; }
      void    padding_right( integer padding ) { m_padding_right = padding; }

      Alignment alignment() const { return m_Align; }
      void      alignment( Alignment align ) { m_Align = align; }

      integer width() const { return m_Width; }
      void    width( integer width ) { m_Width = width; }
    };

    // Dichiarazione anticipata della classe Table
    class Table;

    //!
    //! \brief Represents a cell in a table with alignment, content, and
    //! optional column span.
    //!
    //! The `Cell` class manages the content of a single cell in the table. It
    //! allows specification of alignment, column span, and other properties.
    //!
    class Cell
    {
    private:
      Table *   m_Table    = nullptr;
      string    m_Value    = "";
      Alignment m_Align    = Alignment::LEFT;
      integer   m_col_span = 1;
      integer   m_Width    = 0;

    public:
      //!
      //! \brief Default constructor for an empty cell.
      //!
      Cell() = default;

      //!
      //! \brief Constructs a cell with a value and optional column span.
      //!
      //! \param table Pointer to the table containing the cell.
      //! \param val The string value to be displayed in the cell.
      //! \param col_span The number of columns the cell should span.
      //!
      explicit Cell( Table * table, string_view val = "", integer col_span = 1 );

      string_view value() const { return m_Value; }
      void        value( string_view val ) { m_Value = val; }

      Alignment alignment() const { return m_Align; }
      void      alignment( Alignment const & align ) { m_Align = align; }

      integer col_span() const { return m_col_span; }
      void    col_span( integer col_span ) { m_col_span = col_span; }

      integer width( integer col ) const;
      integer height() const;

      integer maximum_line_width() const;

      string      line( integer idx ) const;
      static void trim_line( std::string & line );

      string render( integer line, integer col ) const;
    };

    //!
    //! \brief Represents a row in a table consisting of multiple cells.
    //!
    //! The `Row` class manages a collection of cells that form a single row in
    //! the table. Each cell in the row can be accessed, modified, and rendered
    //! individually.
    //!
    class Row
    {
    protected:
      using vecCell = vector<Cell>;
      using vecstr  = vector<string>;

      Table * m_Table = nullptr;
      vecCell m_Cells;

    public:
      //!
      //! \brief Default constructor for an empty row.
      //!
      Row() = default;

      //!
      //! \brief Constructs a row with a set of initial cell values.
      //!
      //! \param table Pointer to the table containing the row.
      //! \param cells A vector of strings representing initial cell values.
      //!
      explicit Row( Table * table, vecstr const & cells = vecstr() );

      Table const * table() const { return m_Table; }

      void cells( vecstr const & cells );

      integer num_cells() const { return integer( m_Cells.size() ); }
      integer cell_width( integer idx ) const;
      void    cell_col_span( integer idx, integer span );

      void cell( string_view value );

      Cell const & operator[]( integer idx ) const { return m_Cells[idx]; }
      Cell &       operator[]( integer idx ) { return m_Cells[idx]; }

      integer height() const;

      string render() const;
    };

    //!
    //! \brief The main class for creating and managing a table.
    //!
    //! The `Table` class represents a 2D table structure that supports rows,
    //! cells, and table styles. It provides methods for rendering, alignment,
    //! and other table-related operations.
    //!
    class Table
    {
    public:
      using vecRow    = std::vector<Row>;
      using vecCell   = std::vector<Cell>;
      using vecstr    = std::vector<string>;
      using vecvecstr = std::vector<vecstr>;
      using integer   = int;

    private:
      Style  m_Style;
      string m_Title;
      Row    m_Headings;
      vecRow m_Rows;

      // Helper function to get junction character based on position
      string get_junction_char( integer row_type, integer col, integer nc ) const
      {
        Junction j;

        // Determine connections based on position
        switch ( row_type )
        {
          case 0:  // Top border
            j.up    = false;
            j.down  = true;
            j.left  = ( col == 0 );
            j.right = ( col < nc - 1 );
            break;

          case 1:  // Middle separator (between rows)
            j.up    = true;
            j.down  = true;
            j.left  = ( col == 0 );
            j.right = ( col < nc - 1 );
            break;

          case 2:  // Bottom border
            j.up    = true;
            j.down  = false;
            j.left  = ( col == 0 );
            j.right = ( col < nc - 1 );
            break;

          default:
            j.up    = true;
            j.down  = true;
            j.left  = ( col == 0 );
            j.right = ( col < nc - 1 );
        }

        return resolve_junction( j );
      }

    public:
      //!
      //! \brief Default constructor for an empty table.
      //!
      Table() = default;

      //!
      //! \brief Constructs a table with a given style and initial rows.
      //!
      //! \param style The style object to customize the table's borders and
      //! alignment.
      //! \param rows A 2D vector of strings representing the table's content.
      //!
      explicit Table( Style const & style, vecvecstr const & rows = vecvecstr() ) : m_Style( style )
      {
        this->rows( rows );
      }

      void setup( Style const & style, vecvecstr const & rows = vecvecstr() )
      {
        m_Style = style;
        this->rows( rows );
      }

      void align_column( integer n, Alignment align );
      void add_row( vecstr const & row );

      integer cell_spacing() const;
      integer cell_padding() const;

      vecCell column( integer n ) const;
      integer column_width( integer n ) const;
      integer num_columns() const;

      Style const & style() const { return m_Style; }

      void style( Style const & style ) { m_Style = style; }

      string_view title() const { return m_Title; }

      void title( string_view title ) { m_Title = title; }

      Row const & headings() const { return m_Headings; }

      void headings( vecstr const & headings );

      Row &       row( integer n );
      Row const & row( integer n ) const;

      Row &       operator[]( integer n ) { return this->row( n ); }
      Row const & operator[]( integer n ) const { return this->row( n ); }

      Cell &       operator()( integer i, integer j ) { return ( *this )[i][j]; }
      Cell const & operator()( integer i, integer j ) const { return ( *this )[i][j]; }

      vecRow const & rows() const { return m_Rows; }
      void           rows( vecvecstr const & rows );

      std::string render_separator( integer row_type ) const;

      std::string render() const;
    };

    // =========================================================================
    // IMPLEMENTAZIONI DELLE FUNZIONI DI Cell
    // =========================================================================

    inline Cell::Cell( Table * table, string_view val, integer col_span )
      : m_Table( table ), m_Value( val ), m_col_span( col_span )
    {
      m_Width = utf8_display_width( string( val ) );
    }

    inline integer Cell::width( integer col ) const
    {
      integer const padding{ ( m_col_span - 1 ) * m_Table->cell_spacing() };
      integer       innerWidth{ 0 };

      for ( integer i{ 0 }; i < m_col_span; ++i ) innerWidth += m_Table->column_width( col + i );

      return innerWidth + padding;
    }

    inline integer Cell::height() const
    {
      return static_cast<integer>( std::count( m_Value.begin(), m_Value.end(), '\n' ) + 1 );
    }

    inline integer Cell::maximum_line_width() const
    {
      integer            maxlen{ 0 };
      string             line;
      std::istringstream stream( m_Value );
      while ( std::getline( stream, line ) )
      {
        integer len = utf8_display_width( line );
        if ( len > maxlen ) maxlen = len;
      }
      return maxlen;
    }

    inline string Cell::line( integer idx ) const
    {
      if ( idx < this->height() )
      {
        std::istringstream stream( m_Value );
        string             line;
        for ( integer i{ 0 }; i <= idx; ++i ) std::getline( stream, line );
        trim_line( line );
        return line;
      }
      return "";
    }

    inline void Cell::trim_line( std::string & line )
    {
      auto fun = []( char const c ) -> bool { return std::isspace( static_cast<int>( c ) ) == 0; };
      line.erase( line.begin(), std::find_if( line.begin(), line.end(), fun ) );
      line.erase( std::find_if( line.rbegin(), line.rend(), fun ).base(), line.end() );
    }

    inline string Cell::render( integer line, integer col ) const
    {
      std::stringstream ss;
      integer const     cell_width   = this->width( col );
      string            current_line = this->line( line );
      integer const     line_width   = utf8_display_width( current_line );

      integer const pL{ m_Table->style().padding_left() };
      integer const pR{ m_Table->style().padding_right() };

      integer available_width = cell_width;

      switch ( m_Align )
      {
        case Alignment::LEFT:
        {
          ss << string( pL, ' ' ) << current_line;
          integer remaining = available_width - line_width;
          if ( remaining > 0 )
            ss << string( remaining + pR, ' ' );
          else
            ss << string( pR, ' ' );
          break;
        }
        case Alignment::RIGHT:
        {
          integer remaining = available_width - line_width;
          if ( remaining > 0 )
            ss << string( pL + remaining, ' ' ) << current_line;
          else
            ss << string( pL, ' ' ) << current_line;
          ss << string( pR, ' ' );
          break;
        }
        case Alignment::CENTER:
        {
          integer remaining = available_width - line_width;
          integer left_pad  = remaining / 2;
          integer right_pad = remaining - left_pad;
          ss << string( pL + left_pad, ' ' ) << current_line << string( right_pad + pR, ' ' );
          break;
        }
      }
      return ss.str();
    }

    // =========================================================================
    // IMPLEMENTAZIONI DELLE FUNZIONI DI Row
    // =========================================================================

    inline Row::Row( Table * table, vecstr const & cells ) : m_Table( table )
    {
      for ( auto const & cell : cells ) this->cell( cell );
    }

    inline void Row::cells( vecstr const & cells )
    {
      m_Cells.clear();
      for ( auto const & cell : cells ) this->cell( cell );
    }

    inline integer Row::cell_width( integer idx ) const
    {
      if ( idx < this->num_cells() ) return m_Cells[idx].maximum_line_width();
      return 0;
    }

    inline void Row::cell_col_span( integer idx, integer span )
    {
      if ( span > 0 && idx < this->num_cells() ) m_Cells[idx].col_span( span );
    }

    inline void Row::cell( string_view value )
    {
      m_Cells.emplace_back( m_Table, value );
    }

    inline integer Row::height() const
    {
      integer maxlen = 1;
      for ( auto const & cell : m_Cells )
      {
        if ( cell.height() > maxlen ) maxlen = cell.height();
      }
      return maxlen;
    }

    inline string Row::render() const
    {
      integer const     num_columns{ m_Table->num_columns() };
      integer const     numLines{ this->height() };
      std::stringstream ss;

      Style const   style{ m_Table->style() };
      integer const nc{ static_cast<integer>( m_Cells.size() ) };

      for ( integer l{ 0 }; l < numLines; ++l )
      {
        ss << style.border_left();
        integer c{ 0 };
        while ( c < num_columns )
        {
          if ( c < nc )
          {
            Cell const & cell = m_Cells[c];
            ss << cell.render( l, c );
            if ( cell.col_span() > 1 ) c += cell.col_span() - 1;
          }
          else
          {
            integer col_width = m_Table->column_width( c );
            integer padding   = style.padding_left() + style.padding_right();
            ss << string( col_width + padding, ' ' );
          }

          if ( c < num_columns - 1 ) { ss << style.border_middle(); }
          ++c;
        }
        ss << style.border_right() << '\n';
      }
      return ss.str();
    }

    // =========================================================================
    // IMPLEMENTAZIONI DELLE FUNZIONI DI Table
    // =========================================================================

    inline void Table::align_column( integer n, Alignment align )
    {
      if ( n >= this->num_columns() )
        throw std::out_of_range(
          "Table error: The table just has " + std::to_string( this->num_columns() ) + " columns." );

      // Apply to header if exists
      if ( n < m_Headings.num_cells() ) { m_Headings[n].alignment( align ); }

      // Apply to all rows
      for ( auto & row : m_Rows )
      {
        if ( n < row.num_cells() ) row[n].alignment( align );
      }
    }

    inline void Table::add_row( vecstr const & row )
    {
      m_Rows.emplace_back( this, row );
    }

    inline integer Table::cell_spacing() const
    {
      return this->cell_padding() + 1;
    }

    inline integer Table::cell_padding() const
    {
      return m_Style.padding_left() + m_Style.padding_right();
    }

    inline Table::vecCell Table::column( integer n ) const
    {
      if ( n >= this->num_columns() )
        throw std::out_of_range(
          "Table error: The table just has " + std::to_string( this->num_columns() ) + " columns." );

      vecCell column;
      column.reserve( m_Rows.size() );

      for ( auto const & row : m_Rows )
      {
        if ( n < row.num_cells() )
          column.push_back( row[n] );
        else
          column.push_back( Cell( nullptr, "" ) );
      }

      return column;
    }

    inline integer Table::column_width( integer n ) const
    {
      if ( n >= this->num_columns() )
        throw std::out_of_range(
          "Table error: The table just has " + std::to_string( this->num_columns() ) + " columns." );

      integer maxlen = 0;

      // Check header
      if ( n < m_Headings.num_cells() )
      {
        integer w = m_Headings.cell_width( n );
        if ( w > maxlen ) maxlen = w;
      }

      // Check rows
      for ( auto const & row : m_Rows )
      {
        if ( n < row.num_cells() )
        {
          Cell const & cell = row[n];
          for ( integer i = 0; i < n; ++i )
          {
            if ( i < row.num_cells() )
            {
              Cell const & c = row[i];
              if ( c.col_span() > 1 ) i += c.col_span() - 1;
            }
          }

          integer w = cell.maximum_line_width();
          if ( w > maxlen ) maxlen = w;
        }
      }

      return maxlen;
    }

    inline integer Table::num_columns() const
    {
      integer maxlen = m_Headings.num_cells();
      for ( auto const & row : m_Rows )
      {
        integer cols = 0;
        for ( integer i = 0; i < row.num_cells(); ++i )
        {
          Cell const & cell = row[i];
          cols += cell.col_span();
        }
        if ( cols > maxlen ) maxlen = cols;
      }
      return maxlen;
    }

    inline void Table::headings( vecstr const & headings )
    {
      m_Headings = Row( this, headings );
    }

    inline Row & Table::row( integer n )
    {
      if ( n >= static_cast<integer>( m_Rows.size() ) )
        throw std::out_of_range( "Table error: The table just has " + std::to_string( m_Rows.size() ) + " rows." );
      return m_Rows[n];
    }

    inline Row const & Table::row( integer n ) const
    {
      if ( n >= static_cast<integer>( m_Rows.size() ) )
        throw std::out_of_range( "Table error: The table just has " + std::to_string( m_Rows.size() ) + " rows." );
      return m_Rows[n];
    }

    inline void Table::rows( vecvecstr const & rows )
    {
      m_Rows.clear();
      for ( auto const & row : rows ) m_Rows.emplace_back( this, row );
    }

    inline std::string Table::render_separator( integer row_type ) const
    {
      std::stringstream ss;
      integer const     nc{ this->num_columns() };
      integer const     padding_LR{ m_Style.padding_left() + m_Style.padding_right() };

      for ( integer i{ 0 }; i < nc; ++i )
      {
        string junction_char;
        if ( i == 0 )
        {
          // First column
          if ( row_type == 0 )
            junction_char = m_Style.border_top_left();
          else if ( row_type == 2 )
            junction_char = m_Style.border_bottom_left();
          else
            junction_char = m_Style.border_left_mid();
        }
        else if ( i == nc - 1 )
        {
          // Last column
          if ( row_type == 0 )
            junction_char = m_Style.border_top_right();
          else if ( row_type == 2 )
            junction_char = m_Style.border_bottom_right();
          else
            junction_char = m_Style.border_right_mid();
        }
        else
        {
          // Interior columns
          Junction j;
          j.up          = ( row_type != 0 );
          j.down        = ( row_type != 2 );
          j.left        = true;
          j.right       = true;
          junction_char = resolve_junction( j );
        }

        if ( i == 0 ) ss << junction_char;

        integer const width{ this->column_width( i ) + padding_LR };
        ss << repeat( m_Style.border_mid(), width );

        if ( i < nc - 1 )
        {
          // Junction between columns
          Junction j;
          j.up    = ( row_type != 0 );
          j.down  = ( row_type != 2 );
          j.left  = true;
          j.right = true;
          ss << resolve_junction( j );
        }
        else
        {
          ss << junction_char;
        }
      }
      ss << '\n';
      return ss.str();
    }

    inline std::string Table::render() const
    {
      std::stringstream ss;

      // Calculate total width for title
      integer       innerWidth{ 0 };
      integer const nc{ this->num_columns() };
      integer const padding_LR{ m_Style.padding_left() + m_Style.padding_right() };

      for ( integer c{ 0 }; c < nc; ++c )
      {
        innerWidth += this->column_width( c ) + padding_LR;
        if ( c < nc - 1 ) innerWidth += 1;  // For vertical borders
      }

      // Title
      if ( !m_Title.empty() )
      {
        // Top border
        ss << render_separator( 0 );

        // Title row
        integer title_len  = utf8_display_width( m_Title );
        integer spaceLeft  = ( innerWidth - title_len ) / 2;
        integer spaceRight = innerWidth - title_len - spaceLeft;

        ss << m_Style.border_left() << string( spaceLeft, ' ' ) << m_Title << string( spaceRight, ' ' )
           << m_Style.border_right() << '\n';
      }
      else
      {
        // Top border without title
        ss << render_separator( 0 );
      }

      // Headings
      if ( m_Headings.num_cells() > 0 )
      {
        ss << m_Headings.render();
        ss << render_separator( 1 );
      }

      // Rows
      for ( size_t i = 0; i < m_Rows.size(); ++i )
      {
        ss << m_Rows[i].render();
        if ( i < m_Rows.size() - 1 ) { ss << render_separator( 1 ); }
      }

      // Bottom border
      ss << render_separator( 2 );

      return ss.str();
    }
  }  // namespace Table
}  // namespace Utils

//!
//! \brief Stream insertion operator for rendering a table row to an output
//! stream.
//!
inline Utils::ostream_type & operator<<( Utils::ostream_type & stream, Utils::Table::Row const & row )
{
  return stream << row.render();
}

//!
//! \brief Stream insertion operator for rendering a table to an output stream.
//!
inline Utils::ostream_type & operator<<( Utils::ostream_type & stream, Utils::Table::Table const & table )
{
  return stream << table.render();
}

#endif

//
// eof: Utils_Table.hh
//
