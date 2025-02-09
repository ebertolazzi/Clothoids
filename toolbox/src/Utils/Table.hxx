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
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

//
// eof: Table.hxx
//

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

namespace Utils {

  using std::string;
  using std::string_view;
  using std::vector;

  namespace Table {

    class Table;

    using integer = int;

    //!
    //! \brief Enum class defining alignment types for table cells.
    //!
    //! Provides options for aligning content inside table cells to the LEFT,
    //! RIGHT, or CENTER.
    //!
    using Alignment = enum class Table_align : integer { LEFT, RIGHT, CENTER };

    //!
    //! \brief Defines the style and structure of table borders, padding, and alignment.
    //!
    //! The `Style` class configures visual properties of the table, including border
    //! characters (for all edges and dividers), cell padding, and alignment of text.
    //!
    class Style {
    private:

      char m_border_top          = '-';
      char m_border_top_mid      = '+';
      char m_border_top_left     = '+';
      char m_border_top_right    = '+';

      char m_border_bottom       = '-';
      char m_border_bottom_mid   = '+';
      char m_border_bottom_left  = '+';
      char m_border_bottom_right = '+';

      char m_border_left         = '|';
      char m_border_left_mid     = '+';

      char m_border_mid          = '-';
      char m_border_mid_mid      = '+';

      char m_border_right        = '|';
      char m_border_right_mid    = '+';

      char m_border_middle       = '|';

      integer m_padding_left     = 1;
      integer m_padding_right    = 1;

      Alignment m_Align = Alignment::LEFT;

      integer m_Width = 0;

    public:

      //!
      //! \brief Default constructor initializing the table style with default borders.
      //!
      Style() = default;

      char border_top() const { return m_border_top; }
      void border_top( char borderStyle ) { m_border_top = borderStyle; }

      char border_top_mid() const { return m_border_top_mid; }
      void border_top_mid( char borderStyle ) { m_border_top_mid = borderStyle; }

      char border_top_left() const { return m_border_top_left; }
      void border_top_left( char borderStyle ) { m_border_top_left = borderStyle; }

      char border_top_right() const { return m_border_top_right; }
      void border_top_right( char borderStyle ) { m_border_top_right = borderStyle; }

      char border_bottom() const { return m_border_bottom; }
      void border_bottom( char borderStyle ) { m_border_bottom = borderStyle; }

      char border_bottom_mid() const { return m_border_bottom_mid; }
      void border_bottom_mid( char borderStyle ) { m_border_bottom_mid = borderStyle; }

      char border_bottom_left() const { return m_border_bottom_left; }
      void border_bottom_left( char borderStyle ) { m_border_bottom_left = borderStyle; }

      char border_bottom_right() const { return m_border_bottom_right; }
      void border_bottom_right( char borderStyle) { m_border_bottom_right = borderStyle; }

      char border_left() const { return m_border_left; }
      void border_left( char borderStyle ) { m_border_left = borderStyle; }

      char border_left_mid() const { return m_border_left_mid; }
      void border_left_mid( char borderStyle ) { m_border_left_mid = borderStyle; }

      char border_mid() const { return m_border_mid; }
      void border_mid( char borderStyle ) { m_border_mid = borderStyle; }

      char border_mid_mid() const { return m_border_mid_mid; }
      void border_mid_mid( char borderStyle ) { m_border_mid_mid = borderStyle; }

      char border_right() const { return m_border_right; }
      void border_right( char borderStyle ) { m_border_right = borderStyle; }

      char border_right_mid() const { return m_border_right_mid; }
      void border_right_mid( char borderStyle ) { m_border_right_mid = borderStyle; }

      char border_middle() const { return m_border_middle; }
      void border_middle( char borderStyle ) { m_border_middle = borderStyle; }

      integer padding_left() const { return m_padding_left; }
      void    padding_left( integer padding ) { m_padding_left = padding; }

      integer padding_right() const { return m_padding_right; }
      void    padding_right( integer padding ) { m_padding_right = padding; }

      Alignment alignment() const { return m_Align; }
      void      alignment( Alignment align ) { m_Align = align; }

      integer width() const { return m_Width; }
      void    width( integer width ) { m_Width = width; }
    };

    //!
    //! \brief Represents a cell in a table with alignment, content, and optional column span.
    //!
    //! The `Cell` class manages the content of a single cell in the table. It allows
    //! specification of alignment, column span, and other properties.
    //!
    class Cell {
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
      explicit
      Cell(
        Table*      table,
        string_view val      = "",
        integer     col_span = 1
      );

      string_view value() const { return m_Value; }
      void value( string_view val ) { m_Value = val; }

      Alignment alignment() const { return m_Align; }
      void alignment( Alignment align ) { m_Align = align; }

      integer col_span() const { return m_col_span; }
      void col_span( integer col_span ) { m_col_span = col_span; }

      integer width( integer col ) const;
      integer height() const;

      integer maximum_line_width() const;

      string line( integer idx ) const;
      void trim_line( std::string & line ) const;

      string render( integer line, integer col ) const;
    };

    //!
    //! \brief Represents a row in a table consisting of multiple cells.
    //!
    //! The `Row` class manages a collection of cells that form a single row in the table.
    //! Each cell in the row can be accessed, modified, and rendered individually.
    //!
    class Row {
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
      explicit
      Row(
        Table *        table,
        vecstr const & cells = vecstr()
      );

      Table const * table() const { return m_Table; }

      //vecCell & cells() { return m_Cells; }
      void cells( vecstr const & cells );

      integer num_cells() const { return integer(m_Cells.size()); }
      integer cell_width( integer idx ) const;
      void cell_col_span( integer idx, integer span );

      void cell( string_view value );
      //Cell& cell( integer idx ) { return m_Cells[idx]; }

      Cell const & operator [] ( integer idx ) const { return m_Cells[idx]; }
      Cell       & operator [] ( integer idx )       { return m_Cells[idx]; }

      integer height() const;

      string render() const;
    };

   //!
   //! \brief The main class for creating and managing a table.
   //!
   //! The `Table` class represents a 2D table structure that supports rows, cells,
   //! and table styles. It provides methods for rendering, alignment, and other
   //! table-related operations.
   //!
    class Table {
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

    public:
      //!
      //! \brief Default constructor for an empty table.
      //!
      Table() = default;

      //!
      //! \brief Constructs a table with a given style and initial rows.
      //!
      //! \param style The style object to customize the table's borders and alignment.
      //! \param rows A 2D vector of strings representing the table's content.
      //!
      explicit
      Table(
        Style     const & style,
        vecvecstr const & rows = vecvecstr()
      )
      : m_Style(style) {
        this->rows(rows);
      }

      void
      setup(
        Style     const & style,
        vecvecstr const & rows = vecvecstr()
      ) {
        m_Style = style;
        this->rows(rows);
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

      Row       & row( integer n );
      Row const & row( integer n ) const;

      Row       & operator [] ( integer n )       { return this->row(n); }
      Row const & operator [] ( integer n ) const { return this->row(n); }

      Cell       & operator () ( integer i, integer j )       { return (*this)[i][j]; }
      Cell const & operator () ( integer i, integer j ) const { return (*this)[i][j]; }

      vecRow const & rows() const { return m_Rows; }
      void rows( vecvecstr const & rows );

      std::string
      render_separator(
        char left,
        char mid,
        char right,
        char sep
      ) const;

      std::string render() const;
    };
  }
}

//!
//! \brief Stream insertion operator for rendering a table row to an output stream.
//!
inline
Utils::ostream_type&
operator << ( Utils::ostream_type& stream, Utils::Table::Row const & row ) {
  return stream << row.render();
}

//!
//! \brief Stream insertion operator for rendering a table to an output stream.
//!
inline
Utils::ostream_type&
operator << ( Utils::ostream_type& stream, Utils::Table::Table const & table ) {
  return stream << table.render();
}

//
// eof: Table.hxx
//
