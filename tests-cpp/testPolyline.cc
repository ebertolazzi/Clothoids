//#define _USE_MATH_DEFINES
#include "PolyLine.hh"
#include "Clothoid.hh"
#include <cmath>
#include <iostream>

using G2lib::real_type;

int
main() {

  G2lib::PolyLine      P;
  G2lib::ClothoidCurve C;

  real_type x0  = 0;
  real_type y0  = 0;
  real_type th0 = 0;
  real_type x1  = 0;
  real_type y1  = 3;
  real_type th1 = G2lib::m_pi;

  C.build_G1( x0, y0, th0, x1, y1, th1 );
  C.info(std::cout);

  P.build( C, 0.001 );
  P.info(std::cout);

  std::cout << "\n\nALL DONE FOLKS!!!\n";

  return 0;
}
