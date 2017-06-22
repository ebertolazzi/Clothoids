#define _USE_MATH_DEFINES
#include "Clothoid.hh"
#include <cmath>
#include <iostream>

using Clothoid::valueType ;

int
main() {

  Clothoid::G2solve2arc   g2solve2arc ;
  Clothoid::G2solve3arc   g2solve3arc ;

  valueType m_pi = 3.1415926535897932385 ;

#if 0
  valueType x0  = 0 ;
  valueType y0  = 0 ;
  valueType th0 = m_pi/2+m_pi/4 ;
  valueType k0  = -0.5 ;
  valueType s0  = 0.5 ;

  valueType x1  = 2*sqrt(2.0) ;
  valueType y1  = 2*sqrt(2.0) ;

  valueType th1 = -m_pi/2+m_pi/4 ;
  valueType k1  = -0.5 ;
  valueType s1  = 0.5 ;
#else
  valueType x0  = -1 ;
  valueType y0  = 0 ;
  valueType th0 = m_pi/2 ;
  valueType k0  = -1 ;
  valueType s0  = 0.2 ;

  valueType x1  = 1 ;
  valueType y1  = 0 ;
  valueType th1 = -m_pi/2 ;
  valueType k1  = -1 ;
  valueType s1  = 0.2  ;
#endif

  // test 3 archi
  std::cout << "\n\nTry Solution with THREE arcs\n";
  g2solve3arc.setup( x0, y0, th0, k0, s0, x1, y1, th1, k1, s1 ) ;
  int iter = g2solve3arc.solve() ;
  std::cout << "iter = " << iter << '\n' ;
  
  Clothoid::ClothoidCurve const & S0 = g2solve3arc.getS0() ;
  Clothoid::ClothoidCurve const & S1 = g2solve3arc.getS1() ;
  Clothoid::ClothoidCurve const & SM = g2solve3arc.getSM() ;

  std::cout << "\n\nS0 (NEW)\n" << S0 ;
  std::cout << "\n\nSM (NEW)\n" << SM ;
  std::cout << "\n\nS1 (NEW)\n" << S1 ;

  std::cout << "\n\nALL DONE FOLKS!!!\n" ;

  return 0 ;
}
