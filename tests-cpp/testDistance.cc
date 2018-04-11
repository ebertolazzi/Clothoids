#define _USE_MATH_DEFINES
#include "Clothoid.hh"
#include <cmath>
#include <iostream>

using G2lib::valueType ;

int
main() {

  G2lib::ClothoidCurve C ;
  valueType x0     = 0 ;
  valueType y0     = 2 ;
  valueType theta0 = 0 ;
  valueType kappa0 = 10 ;
  valueType dk     = -1 ;// 0 ;
  valueType L      = 10 ;
  C.build( x0, y0, theta0, kappa0, dk, L ) ;

  //valueType x = -2.75 ;
  //valueType y = -3.2000000000000001776 ;
  valueType x = 10 ;
  valueType y = 15 ;
  valueType X, Y, S ;

  valueType d = C.closestPoint( x, y, X, Y, S );

  std::cout << "\nd = " << d
            << "\nX = " << X
            << "\nY = " << Y
            << "\nS = " << S
            << "\n" ;

  std::cout << "\n\nALL DONE FOLKS!!!\n" ;

  return 0 ;
}
