#define _USE_MATH_DEFINES
#include "Clothoid.hh"
#include <cmath>
#include <iostream>

using G2lib::valueType ;

int
main() {

  G2lib::G2solve2arc g2solve2arc ;
  G2lib::G2solve3arc g2solve3arc ;

  valueType m_pi = 3.1415926535897932385 ;

#if 1
  valueType x0  = -1 ;
  valueType y0  =  0 ;
  valueType th0 = -2 ;
  valueType k0  = 0.909297426825682 ;
  valueType x1  = 1 ;
  valueType y1  = 0 ;
  valueType th1 = 2 ;
  valueType k1  = 0.909297426825682 ;
#else
  valueType x0  = -1 ;
  valueType y0  = 0  ;
  valueType th0 = m_pi*0.9;
  valueType k0  = 0.2 + 1e-10 ;

  valueType x1  = 1 ;
  valueType y1  = 0 ;
  valueType th1 = -m_pi ;
  valueType k1  = 0.2 + 0 ;
#endif

  // test 3 archi
  std::cout << "\n\nTry Solution with THREE arcs\n";
  int iter = g2solve3arc.build( x0, y0, th0, k0, x1, y1, th1, k1 ) ;
  std::cout << "iter = " << iter << '\n' ;
  
  G2lib::ClothoidCurve const & S0 = g2solve3arc.getS0() ;
  G2lib::ClothoidCurve const & S1 = g2solve3arc.getS1() ;
  G2lib::ClothoidCurve const & SM = g2solve3arc.getSM() ;

  std::cout << "\n\nS0 (NEW)\n" << S0 ;
  std::cout << "\n\nSM (NEW)\n" << SM ;
  std::cout << "\n\nS1 (NEW)\n" << S1 ;

  std::cout
    << "\nx  = " << S0.xEnd()       << " " << SM.xBegin()     << " err = " << S0.xEnd()-SM.xBegin()
    << "\ny  = " << S0.yEnd()       << " " << SM.yBegin()     << " err = " << S0.yEnd()-SM.yBegin()
    << "\nth = " << S0.thetaEnd()   << " " << SM.thetaBegin() << " err = " << S0.thetaEnd()-SM.thetaBegin()
    << "\nx  = " << S1.xBegin()     << " " << SM.xEnd()       << " err = " << S1.xBegin()-SM.xEnd()
    << "\ny  = " << S1.yBegin()     << " " << SM.yEnd()       << " err = " << S1.yBegin()-SM.yEnd()
    << "\nth = " << S1.thetaBegin() << " " << SM.thetaEnd()   << " err = " << S1.thetaBegin()-SM.thetaEnd()
    << '\n' ;

  std::cout << "\n\nALL DONE FOLKS!!!\n" ;

  return 0 ;
}
