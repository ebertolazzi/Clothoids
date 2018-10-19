//#define _USE_MATH_DEFINES
#include "Triangle2D.hh"
#include <cmath>
#include <iostream>

using G2lib::real_type;
using G2lib::int_type;

int
main() {

  // check constructors
  real_type x0     = 0;
  real_type y0     = 2;
  real_type x1     = 4;
  real_type y1     = 3.5;
  real_type x2     = 6;
  real_type y2     = 1;


  G2lib::Triangle2D T1;
  T1.build( x0, y0, x1, y1, x2, y2 );

  int_type icode = T1.isInside( T1.baricenterX(), T1.baricenterY() );

  std::cout << "icode = " << icode << "\n";

  icode = T1.isInside( T1.P1() );

  std::cout << "icode = " << icode << "\n";

  icode = T1.isInside( 10, 20 );

  std::cout << "icode = " << icode << "\n";



  std::cout << "\n\nALL DONE FOLKS!!!\n";

  return 0;
}
