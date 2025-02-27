#include "Clothoids.hh"
#include "Clothoids_fmt.hh"

namespace G2lib {

  AsyPlot::AsyPlot( string _filename, bool _showAxes )
  : filename(_filename)
  , showAxes(_showAxes)
  {
    UTILS_ASSERT( openFile(), "Failed to open file {}\n", filename );
    initFile();
  }

  AsyPlot::~AsyPlot() {
    if ( showAxes ) displayAxes();
    if ( closeFile() ) compileFile();
  }

  void
  AsyPlot::compileFile() {
    string const cmdComp = "asy -f pdf "+ filename;
    system(cmdComp.c_str());
    string const pdfFile = filename.substr(0,filename.find(".asy")) + ".pdf";
    std::cout << pdfFile << std::endl;
    string const cmdOpen = "(okular " + pdfFile + " &> /dev/null )&";
    system(cmdOpen.c_str());
  }

  void
  AsyPlot::initFile() {
    file
      << "// File generated automatically from C++ \n\n\n"
      << "import graph;\n"
      << "include \"clothoidLib.asylib\";\n"
      << "size(14cm,7cm);\n"
      << "\n\n\n";
  }

  void
  AsyPlot::drawClothoid( ClothoidCurve const & c,
                         string_view   const   penna,
                         real_type     const   offset ) const {
  	if (offset == 0.) {
      fmt::print( file,
        "path pclot = clothoidPoints(({},{}),{},{},{},{},100,0);\n"
        "pen  penna = {};\n"
        "draw(pclot, penna);\n\n",
        c.x_begin(), c.y_begin(),
        c.theta_begin(), c.kappa_begin(), c.dkappa(), c.length()
      );
	  } else {
      fmt::print( file,
        "path pclot = clothoidOffset(({},{}),{},{},{},{},100,{});\n"
        "pen  penna = {};\n"
        "draw(pclot, penna);\n\n",
		    c.x_begin(), c.y_begin(),
        c.theta_begin(), c.kappa_begin(), c.dkappa(), c.length(),
        offset
      );
    }
  }

  void
  AsyPlot::dot( real_type const x, real_type const y, string_view const penna ) const {
    file << "dot((" << x << "," << y << ")," << penna << ");\n\n";
  }

  void
  AsyPlot::triangle( Triangle2D const & t, string_view const penna ) const {
    file
      << "draw((" << t.x1() << "," << t.y1() << ") -- "
      << '(' << t.x2() << ',' << t.y2() << ") -- "
      << '(' << t.x3() << ',' << t.y3() << ") -- cycle, "
      << penna << ");\n\n";
  }

  void
  AsyPlot::drawRect(
    real_type const x0, real_type const y0,
    real_type const x1, real_type const y1,
    real_type const x2, real_type const y2,
    real_type const x3, real_type const y3,
    string_view const penna
  ) const {
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
  AsyPlot::displayAxes(
    string_view const labX,
    string_view const labY,
    real_type   const xmin,
    real_type   const xmax,
    real_type   const ymin,
    real_type   const ymax
  ) const {
	  file
      << "xaxis(\"" << labX << "\", black+fontsize(7pt),xmin="
      << xmin << ",xmax=" << xmax
      << ",Ticks(Step=2,step=1,NoZero,Size=.8mm, size=.4mm));\n"
      << "yaxis(\"" << labY << "\", black+fontsize(7pt),ymin="
      << ymin << ",ymax=" << ymax
      << ",Ticks(Step=2,step=1,NoZero,Size=.8mm, size=.4mm));\n";
  }

  void
  AsyPlot::drawLine(
    real_type const x0, real_type const y0,
    real_type const x1, real_type const y1,
    string_view const penna
  ) const {
    file
      << "draw((" << x0 << "," << y0 << ") -- "
      << "(" << x1 << "," << y1 << "), " << penna << ");\n\n";
  }

  void
  AsyPlot::label(
    string_view const text,
    real_type   const x,
    real_type   const y,
    string_view const placement,
    string_view const penna
  ) const {
  	file << "label(\"" << text << "\", (" << x << ", " << y << "), ";
    if ( !placement.empty() ) file << placement << ", ";
	  file << penna << ");\n\n";
  }

  bool
  AsyPlot::openFile() {
    file.open(filename.c_str());
    return file.is_open();
  }

  bool
  AsyPlot::closeFile() {
    file.close();
    return !file.is_open();
  }

}
