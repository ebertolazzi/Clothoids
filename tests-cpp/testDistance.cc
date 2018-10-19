//#define _USE_MATH_DEFINES
#include "Clothoid.hh"
#include <cmath>
#include <iostream>

using G2lib::real_type;

int
main() {

  G2lib::ClothoidCurve C;
  real_type x0     = 0;
  real_type y0     = 2;
  real_type theta0 = 0;
  real_type kappa0 = 10;
  real_type dk     = -1;// 0;
  real_type L      = 10;
  C.build( x0, y0, theta0, kappa0, dk, L );

  //real_type x = -2.75;
  //real_type y = -3.2000000000000001776;
  real_type x = 10;
  real_type y = 15;
  real_type X, Y, S;

  real_type d = C.closestPoint( x, y, X, Y, S );

  std::cout << "\nd = " << d
            << "\nX = " << X
            << "\nY = " << Y
            << "\nS = " << S
            << "\n";

  std::cout << "\n\nALL DONE FOLKS!!!\n";

  return 0;
}
