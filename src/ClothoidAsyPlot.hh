/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                | 
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi and Marco Frego                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |      email: marco.frego@unitn.it                                         |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#ifndef ASYPLOT_H
#define ASYPLOT_H

#include "Clothoid.hh"

#include <string>
#include <fstream>
#include <iostream>
#include <cmath>

//! Clothoid computations routine
namespace Clothoid {

  using std::string;
  using std::ofstream;

  class AsyPlot {
  public:
    AsyPlot( string filename, bool showAxes );
    ~AsyPlot();

    void
    drawClothoid( Clothoid::ClothoidCurve const& c,
                  std::string const & penna="black",
                  valueType offset = 0.) const ;
    
    void dot( valueType x, valueType y, string const & penna="black" ) const ;
    void triangle(const Triangle2D<valueType> & t, string const & penna="black" ) const ;

    void
    drawRect( valueType x0, valueType y0,
              valueType x1, valueType y1,
              valueType x2, valueType y2,
              valueType x3, valueType y3,
              string const & penna="black") const ;

    void
    drawLine( valueType x0, valueType y0,
              valueType x1, valueType y1,
              std::string const & penna="black" ) const ;
    void
    label( string const & text,
           valueType      x,
           valueType      y,
           string const & placement = "",
           string const & penna = "black" ) const ;

    void
    displayAxes( string const & labX,
                 string const & labY,
                 valueType      xmin,
                 valueType      xmax,
                 valueType      ymin,
                 valueType      ymax ) const ;

  private:
    mutable ofstream file;
    string  filename;
    bool showAxes;
    bool openFile();
    bool closeFile();
    void initFile();
    void displayAxes() const ;
    void compileFile();
  };
}

#endif 
