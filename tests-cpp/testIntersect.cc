#define _USE_MATH_DEFINES
#include "Clothoid.hh"
#include <cmath>
#include <iostream>

using G2lib::valueType ;

int
main() {

  G2lib::ClothoidCurve C0, C1 ;
  valueType x0     = 0 ;
  valueType y0     = 0 ;
  valueType theta0 = G2lib::m_pi*0.7 ;
  valueType k0     = -0.1 ;
  valueType dk0    = 0.05 ;
  valueType L0     = 20 ;
  C0.build( x0, y0, theta0, k0, dk0, L0 );

  valueType x1     = -10 ;
  valueType y1     = 0 ;
  valueType theta1 = G2lib::m_pi/4 ;
  valueType k1     = 0.1 ;
  valueType dk1    = -0.05 ;
  valueType L1     = 20 ;
  C1.build( x1, y1, theta1, k1, dk1, L1 );

  G2lib::vector<G2lib::valueType> s1, s2 ;
  G2lib::indexType max_iter  = 10 ;
  G2lib::valueType tolerance = 1e-8 ;

  C0.intersect( C1, s1, s2, max_iter, tolerance );

  std::cout << "L0 = " << C0.getL() << '\n' ;
  std::cout << "L1 = " << C1.getL() << '\n' ;

  for ( int i = 0 ; i < s1.size() ; ++i )
    std::cout << "s1[ " << i << "] = " << s1[i] << '\n' ;

  for ( int i = 0 ; i < s2.size() ; ++i )
    std::cout << "s2[ " << i << "] = " << s2[i] << '\n' ;

  for ( int i = 0 ; i < s1.size() ; ++i )
    std::cout << "x["  << i << "] = " << C0.X(s1[i])
              << " y[" << i << "] = " << C0.Y(s1[i])
              << '\n' ;

  for ( int i = 0 ; i < s2.size() ; ++i )
    std::cout << "x["  << i << "] = " << C1.X(s2[i])
              << " y[" << i << "] = " << C1.Y(s2[i])
              << '\n' ;

  std::cout << "\n\nALL DONE FOLKS!!!\n" ;

  return 0 ;
}
