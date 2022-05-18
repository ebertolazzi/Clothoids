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
 |      Via Sommarive 9, I-38123 Povo, Trento, Italy                        |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

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


///
/// file: Table.cc
///

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include "Utils.hh"

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cctype>
#include <algorithm>
#include <functional>
#include <stdexcept>

namespace Utils {

  namespace Table {

    using std::istringstream;
    using std::stringstream;
    using std::count;
    using std::getline;
    using std::find_if;
    using std::isspace;
    using std::setw;
    using std::left;
    using std::right;
    using std::setfill;
    using std::bind;
    using std::out_of_range;
    using std::to_string;
    using std::transform;

    Cell::Cell(
      Table*         table,
      string const & val,
      integer        col_span
    )
    : m_Table(table)
    , m_Value(val)
    , m_col_span(col_span)
    {
      m_Width = integer(val.length());
    }

    integer
    Cell::width( integer col ) const {
      integer padding    = (m_col_span - 1) * m_Table->cell_spacing();
      integer innerWidth = 0;

      for( integer i = 0; i < m_col_span; ++i )
        innerWidth += m_Table->column_width(col + i);

      return innerWidth + padding;
    }

    integer
    Cell::height() const {
      return integer( count( m_Value.begin(), m_Value.end(), '\n') + 1);
    }

    integer
    Cell::maximum_line_width() const {
      integer maxlen = 0;
      string  line;
      istringstream stream( m_Value );
      while( getline(stream, line) ) {
        integer len = integer( line.length() );
        if ( len > maxlen ) maxlen = len;
      }
      return maxlen;
    }

    string
    Cell::line( integer idx ) const {
      if ( idx < this->height() ) {
        istringstream stream(m_Value);
        string line;
        for( integer i = 0; i <= idx; ++i) getline(stream, line);
        this->trim_line(line);
        return line;
      }
      return "";
    }

    void
    Cell::trim_line( string & line ) const {
      auto fun = [] ( char c ) -> bool { return isspace(int(c)) == 0; };
      line.erase( line.begin(), find_if( line.begin(), line.end(), fun ) );
      line.erase( find_if( line.rbegin(), line.rend(), fun ).base(), line.end() );
    }

    string
    Cell::render( integer line, integer col ) const {
      stringstream ss;
      integer width = this->width(col);
      integer pL    = m_Table->style().padding_left();
      integer pR    = m_Table->style().padding_right();

      switch( m_Align ) {
      case Alignment::LEFT:
        ss << string( pL, ' ' )
           << setw(width) << left << setfill(' ') << this->line(line)
           << string( pR, ' ' );
        break;
      case Alignment::RIGHT:
        ss << string( pL, ' ' )
           << setw(width) << right << setfill(' ') << this->line(line)
           << string( pR, ' ' );
        break;
      case Alignment::CENTER:
        {
          string  val        = this->line(line);
          integer innerWidth = width + m_Table->cell_padding();
          integer spaceLeft  = (innerWidth-integer(val.length()))/2;
          ss << string(spaceLeft, ' ')
             << setw(innerWidth - spaceLeft) << left << setfill(' ') << this->line(line);
        }
        break;
      }
      return ss.str();
    }

    Row::Row( Table * table, vecstr const & cells )
    : m_Table(table) {
      for_each(
        cells.begin(), cells.end(),
        bind(
          static_cast<void(Row::*)(string const &)>(&Row::cell),
          this,
          std::placeholders::_1
        )
      );
    }

    void
    Row::cells( vecstr const & cells ) {
      for_each(
        cells.begin(), cells.end(),
        bind(
          static_cast<void(Row::*)(string const &)>(&Row::cell),
          this,
          std::placeholders::_1
        )
      );
    }

    integer
    Row::cell_width( integer idx ) const {
      if ( idx < this->num_cells() ) return m_Cells[idx].maximum_line_width();
      return 0;
    }

    void
    Row::cell_col_span( integer idx, integer span ) {
      if ( span > 0 && idx < this->num_cells() ) m_Cells[idx].col_span(span);
    }

    void
    Row::cell( string const & value ) {
      m_Cells.push_back( Cell( m_Table, value ) );
    }

    integer
    Row::height() const {
      integer maxlen = 1;
      for_each(
        m_Cells.begin(), m_Cells.end(),
        [&maxlen]( Cell const & cell ) -> void {
          if ( cell.height() > maxlen ) maxlen = cell.height();
        }
      );
      return maxlen;
    }

    string
    Row::render() const {
      integer num_columns = m_Table->num_columns();
      integer paddingLR   = m_Table->style().padding_left()+
                            m_Table->style().padding_right();
      integer numLines    = this->height();
      stringstream ss;

      Style style = m_Table->style();
      integer nc = integer(m_Cells.size());

      for ( integer l = 0; l < numLines; ++l ) {
        ss << style.border_left();
        for( integer c = 0; c < num_columns; ++c ) {
          if ( c < nc ) {
            Cell const & C = m_Cells[c];
            ss << C.render(l, c);
            if ( C.col_span() > 1 ) c += C.col_span() - 1;
          } else {
            ss << string(m_Table->column_width(c)+paddingLR,' ');
          }
          if ( c < num_columns-1 ) ss << style.border_middle();
        }
        ss << style.border_right() << '\n';
      }
      return ss.str();
    }

    void
    Table::align_column( integer n, Alignment align ) {
      if ( n > this->num_columns() )
        throw out_of_range(
          "Table error: The table just has " + to_string(this->num_columns()) + " columns."
        );
      for_each(
        m_Rows.begin(), m_Rows.end(),
        [n, align]( Row & row ) -> void {
          if ( n < row.num_cells() ) row[n].alignment(align);
        }
      );
    }

    void
    Table::add_row( vecstr const & row ) {
      m_Rows.push_back( Row(this, row) );
    }

    integer
    Table::cell_spacing() const {
      return this->cell_padding() + 1;
    }

    integer
    Table::cell_padding() const {
      return m_Style.padding_left() + m_Style.padding_right();
    }

    typename Table::vecCell
    Table::column( integer n ) const {
      if ( n > this->num_columns() )
        throw out_of_range(
          "Table error: The table just has " + to_string(this->num_columns()) + " columns."
        );

      vecCell column(m_Rows.size());

      transform(
        m_Rows.begin(), m_Rows.end(), column.begin(),
        [n]( Row const & row ) -> Cell { return row[n]; }
      );

      return column;
    }

    integer
    Table::column_width( integer n ) const {
      if ( n > this->num_columns() )
        throw out_of_range(
          "Table error: The table just has " + to_string(this->num_columns()) + " columns."
        );
      integer maxlen = 0;
      auto fun = [&maxlen, n]( Row const & row ) -> void {
        if ( n < row.num_cells() ) {
          integer cw = row.cell_width(n);
          if ( cw > maxlen ) maxlen = cw;
        }
      };
      fun( m_Headings );
      for_each( m_Rows.begin(), m_Rows.end(), fun );
      return maxlen;
    }

    integer
    Table::num_columns() const {
      integer maxlen = m_Headings.num_cells();
      auto fun = [&maxlen]( Row const & row ) -> void {
        integer nc = row.num_cells();
        if ( nc > maxlen ) maxlen = nc;
      };
      for_each( m_Rows.begin(), m_Rows.end(), fun );
      return maxlen;
    }

    void
    Table::headings( vecstr const & headings ) {
      m_Headings = Row(this, headings);
    }

    Row &
    Table::row( integer n ) {
      if ( n >= integer(m_Rows.size()) )
        throw out_of_range(
          "Table error: The table just has " + to_string(m_Rows.size()) + " rows."
        );
      return m_Rows[n];
    }

    Row const &
    Table::row( integer n ) const {
      if ( n >= integer(m_Rows.size()) )
        throw out_of_range(
          "Table error: The table just has " + to_string(m_Rows.size()) + " rows."
        );
      return m_Rows[n];
    }

    void
    Table::rows( vecvecstr const & rows ) {
      m_Rows = vecRow();
      for_each(
        rows.begin(), rows.end(),
        [this]( vecstr const & row ) {
          m_Rows.push_back( Row(this, row) );
        }
      );
    }

    string
    Table::render_separator(
      char left,
      char mid,
      char right,
      char sep
    ) const {
      stringstream ss;
      ss << left;
      integer padding_LR = m_Style.padding_left() + m_Style.padding_right();
      integer nc         = this->num_columns();
      for ( integer i = 0; i < nc; ++i ) {
        integer width = this->column_width(i) + padding_LR;
        for ( integer j = 0; j < width; ++j ) ss << sep;
        if ( i+1 < nc ) ss << mid;
        else            ss << right;
      }
      ss << '\n';
      return ss.str();
    }

    string
    Table::render() const {
      stringstream ss;
      string sep = this->render_separator(
        m_Style.border_left_mid(),
        m_Style.border_mid_mid(),
        m_Style.border_right_mid(),
        m_Style.border_mid()
      );

      if ( m_Title.length() > 0 ) {
        integer innerWidth = (this->num_columns() - 1) * this->cell_spacing() + this->cell_padding();
        for ( integer c = 0; c < this->num_columns(); ++c )
          innerWidth += this->column_width(c);

        integer spaceLeft = (innerWidth - integer(m_Title.length())) / 2;

        ss << m_Style.border_top_left()
           << string(innerWidth, m_Style.border_top())
           << m_Style.border_top_right()
           << '\n'
           << m_Style.border_left()
           << string(spaceLeft, ' ')
           << left << setw(innerWidth - spaceLeft) << setfill(' ') << m_Title
           << m_Style.border_right()
           << '\n' << sep;
      } else {
        ss << render_separator(
          m_Style.border_top_left(),
          m_Style.border_top_mid(),
          m_Style.border_top_right(),
          m_Style.border_top()
        );
      }

      if ( m_Headings.num_cells() > 0 )
        ss << m_Headings.render() << sep;

      if ( !m_Rows.empty() ) {
        for_each(
          m_Rows.begin(), --m_Rows.end(),
          [&ss, sep]( Row const & row ) -> void {
            if ( row.num_cells() > 0 ) ss << row.render() << sep;
          }
        );
        ss << m_Rows.back().render();
      }

      ss << this->render_separator(
              m_Style.border_bottom_left(),
              m_Style.border_bottom_mid(),
              m_Style.border_bottom_right(),
              m_Style.border_bottom()
            );
      return ss.str();
    }
  }
}

#endif

///
/// eof: Table.cc
///
