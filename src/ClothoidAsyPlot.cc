#include "Clothoid.hh"
#include "Triangle2D.hh"
#include "ClothoidAsyPlot.hh"

namespace Clothoid {

  AsyPlot::AsyPlot( string _filename, bool _showAxes )
  : filename(_filename)
  , showAxes(_showAxes)
  {
    G2LIB_ASSERT( openFile(), "Failed to open file " << filename )
    initFile();
  }

  AsyPlot::~AsyPlot() {
    if ( showAxes ) displayAxes();
    if ( closeFile() ) compileFile();
  }

  void
  AsyPlot::compileFile() {
    string cmdComp = "asy -f pdf "+ filename;
    system(cmdComp.c_str());
    string pdfFile = filename.substr(0,filename.find(".asy")) + ".pdf";
    std::cout << pdfFile << std::endl;
    string cmdOpen = "(okular " + pdfFile + " &> /dev/null )&";
    system(cmdOpen.c_str());
  }

  void
  AsyPlot::initFile() {
    file
      << "// File generated automatically from C++ \n\n\n"
      << "import graph;\n"
      << "include \"clothoidLib.asylib\";\n"
      << "size(14cm,7cm);\n"
      << "\n\n\n" ;
  }

  void
  AsyPlot::drawClothoid( ClothoidCurve const & c,
                         string const & penna,
                         valueType offset ) const {
  	if (offset == 0.) {
      file
        << "path pclot = clothoidPoints(("
        << c.getX0()      << ','
        << c.getY0()      << "),"
        << c.getTheta0()  << ','
        << c.getKappa()   << ','
        << c.getKappa_D() << ','
        << c.getL()       << ','
        << "100,0);\n"
        << "pen penna = " << penna << ";\n"
        << "draw(pclot, penna);\n\n";
	  } else {
      file
        << "path pclot = clothoidOffset(("
		    << c.getX0() << ','
        << c.getY0() << "),"
        << c.getTheta0() << ','
	      << c.getKappa() << ','
        << c.getKappa_D() << ','
        << c.getL() << ','
        << "100,"
        << offset
        << "); \n"
        << "pen penna = " << penna << ";\n"
        << "draw(pclot, penna);\n\n";
    }
  }

  void
  AsyPlot::dot( valueType x, valueType y, string const & penna ) const {
    file << "dot((" << x << "," << y << ")," << penna << ");\n\n";
  }

  void
  AsyPlot::triangle( T2D const & t, string const & penna ) const {
    file
      << "draw((" << t.x1() << "," << t.y1() << ") -- "
      << '(' << t.x2() << ',' << t.y2() << ") -- "
      << '(' << t.x3() << ',' << t.y3() << ") -- cycle, "
      << penna << ");\n\n";
  }

  void
  AsyPlot::drawRect( valueType x0, valueType y0,
                     valueType x1, valueType y1,
                     valueType x2, valueType y2,
                     valueType x3, valueType y3,
                     string const & penna ) const {
	file
    << "fill((" << x0 << "," << y0 << ") -- "
    << '(' << x1 << ',' << y1 << ") -- "
    << '(' << x2 << ',' << y2 << ") -- "
    << '(' << x3 << ',' << y3 << ") -- cycle, "
    << penna << ");\n\n";
  }

  void
  AsyPlot::displayAxes() const {
    file
      << "xaxis(\"$x$\", black+fontsize(7pt),xmin=0,xmax=5,Ticks(Step=2,step=1,NoZero,Size=.8mm, size=.4mm));\n"
      << "yaxis(\"$y$\", black+fontsize(7pt),ymin=0,ymax=5,Ticks(Step=2,step=1,NoZero,Size=.8mm, size=.4mm));\n";
  }

  void
  AsyPlot::displayAxes( string const & labX,
                        string const & labY,
                        valueType xmin,
                        valueType xmax,
                        valueType ymin,
                        valueType ymax ) const {
	  file
      << "xaxis(\"" << labX << "\", black+fontsize(7pt),xmin="
      << xmin << ",xmax=" << xmax
      << ",Ticks(Step=2,step=1,NoZero,Size=.8mm, size=.4mm));\n"
      << "yaxis(\"" << labY << "\", black+fontsize(7pt),ymin="
      << ymin << ",ymax=" << ymax
      << ",Ticks(Step=2,step=1,NoZero,Size=.8mm, size=.4mm));\n";
  }

  void
  AsyPlot::drawLine( valueType x0, valueType y0,
                     valueType x1, valueType y1,
                     string const & penna ) const {
    file
      << "draw((" << x0 << "," << y0 << ") -- "
      << "(" << x1 << "," << y1 << "), " << penna << ");\n\n";
  }

  void
  AsyPlot::label( string const & text,
                  valueType      x,
                  valueType      y,
                  string const & placement,
                  string const & penna ) const {
  	file
      << "label(\"" << text << "\", (" << x << ", " << y << "), "
  	  << (!placement.empty() ? placement + ", " : "")
	    << penna << ");\n\n";
  }

  bool
  AsyPlot::openFile() {
    file.open(filename);
    return file.is_open();
  }

  bool
  AsyPlot::closeFile() {
    file.close();
    return !file.is_open();
  }

}
