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
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |      email: marco.frego@unitn.it                                         |
 |                                                                          |
\*--------------------------------------------------------------------------*/

namespace G2lib {

  using std::string;
  using std::ofstream;

  class AsyPlot {
  public:
    AsyPlot( string filename, bool showAxes );
    ~AsyPlot();

    void
    drawClothoid(
      ClothoidCurve const & c,
      string_view           penna  = "black",
      real_type             offset = 0
    ) const;

    void dot( real_type x, real_type y, string_view penna="black" ) const;
    void triangle(Triangle2D const & t, string_view penna="black" ) const;

    void
    drawRect(
      real_type x0, real_type y0,
      real_type x1, real_type y1,
      real_type x2, real_type y2,
      real_type x3, real_type y3,
      string_view penna = "black"
    ) const;

    void
    drawLine(
      real_type x0, real_type y0,
      real_type x1, real_type y1,
      string_view penna = "black"
    ) const;

    void
    label(
      string_view text,
      real_type   x,
      real_type   y,
      string_view placement = "",
      string_view penna     = "black"
    ) const;

    void
    displayAxes(
      string_view labX,
      string_view labY,
      real_type   xmin,
      real_type   xmax,
      real_type   ymin,
      real_type   ymax
    ) const;

  private:
    mutable ofstream file;
    string  filename;
    bool    showAxes;
    bool    openFile();
    bool    closeFile();
    void    initFile();
    void    displayAxes() const;
    void    compileFile();
  };
}
