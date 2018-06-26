#define _USE_MATH_DEFINES
#include "Clothoid.hh"
#include <cmath>
#include <iostream>

using G2lib::real_type ;

int
main() {

  G2lib::ClothoidCurve C0, C1 ;
  real_type x0     = 0 ;
  real_type y0     = 0 ;
  real_type theta0 = G2lib::m_pi*0.7 ;
  real_type k0     = -0.1 ;
  real_type dk0    = 0.05 ;
  real_type L0     = 20 ;
  C0.build( x0, y0, theta0, k0, dk0, L0 );

  real_type x1     = -10 ;
  real_type y1     = 0 ;
  real_type theta1 = G2lib::m_pi/4 ;
  real_type k1     = 0.1 ;
  real_type dk1    = -0.05 ;
  real_type L1     = 20 ;
  C1.build( x1, y1, theta1, k1, dk1, L1 );

  std::vector<G2lib::real_type> s1, s2 ;
  G2lib::int_type  max_iter  = 10 ;
  G2lib::real_type tolerance = 1e-8 ;

  C0.intersect( C1, s1, s2, max_iter, tolerance );

  std::cout << "L0 = " << C0.length() << '\n' ;
  std::cout << "L1 = " << C1.length() << '\n' ;

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
