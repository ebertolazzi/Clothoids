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

  valueType x0  = 0 ;
  valueType y0  = 0  ;
  valueType th0 = 2.9095201640181596 ;
  valueType k0  = 0.10000000000000001 ;

  valueType x1  = 1 ;
  valueType y1  = 0 ;
  valueType th1 = -3.1101767270538954 ;
  valueType k1  = 1 ;

  // test 3 archi
  std::cout << "\n\nTry Solution with THREE arcs\n";
  int iter = g2solve3arc.build( x0, y0, th0, k0, x1, y1, th1, k1 ) ;
  std::cout << "iter = " << iter << '\n' ;
  
  Clothoid::ClothoidCurve const & S0 = g2solve3arc.getS0() ;
  Clothoid::ClothoidCurve const & S1 = g2solve3arc.getS1() ;
  Clothoid::ClothoidCurve const & SM = g2solve3arc.getSM() ;

  std::cout << "\n\nS0 (NEW)\n" << S0 ;
  std::cout << "\n\nSM (NEW)\n" << SM ;
  std::cout << "\n\nS1 (NEW)\n" << S1 ;

  std::cout
    << "\nx  = " << S0.Xend()       << " " << SM.Xbegin()     << " err = " << S0.Xend()-SM.Xbegin()
    << "\ny  = " << S0.Yend()       << " " << SM.Ybegin()     << " err = " << S0.Yend()-SM.Ybegin()
    << "\nth = " << S0.ThetaEnd()   << " " << SM.ThetaBegin() << " err = " << S0.ThetaEnd()-SM.ThetaBegin()
    << "\nx  = " << S1.Xbegin()     << " " << SM.Xend()       << " err = " << S1.Xbegin()-SM.Xend()
    << "\ny  = " << S1.Ybegin()     << " " << SM.Yend()       << " err = " << S1.Ybegin()-SM.Yend()
    << "\nth = " << S1.ThetaBegin() << " " << SM.ThetaEnd()   << " err = " << S1.ThetaBegin()-SM.ThetaEnd()
    << '\n' ;

  std::cout << "\n\nALL DONE FOLKS!!!\n" ;

  return 0 ;
}
