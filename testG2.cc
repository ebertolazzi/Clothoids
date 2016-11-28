#define _USE_MATH_DEFINES
#include "Clothoid.hh"
#include <cmath>
#include <iostream>

using Clothoid::valueType ;

int
main() {

  Clothoid::G2solve2arc   g2solve2arc ;
  Clothoid::G2solve3arc   g2solve3arc ;

  valueType x0  = 0 ;
  valueType y0  = 0 ;
  valueType th0 = 0 ;
  valueType k0  = -0.2 ;
  valueType f0A = 0.33 ;
  valueType x1  = 1 ;
  valueType y1  = 0 ;
  valueType th1 = 0 ;
  valueType k1  = 0.5 ;
  valueType f1B = 0.33 ;
  
  x0 = -1 ;
  y0 = 0 ;
  th0 = M_PI/3 ;
  k0  = 0.1*1.2 ;
  x1 = 1 ;
  y1 = 0 ;
  th1 = M_PI/10 ;
  k1 = -0.1*0.5 ;
  f0A = 0.3 ;
  f1B = 0.4 ;
  /*
  x0  = -2 ;
  y0  = 3 ;
  th0 = 1.047197551196598 ;
  k0  = 10.199999999999999 ;
  f0A = 0.322000000000000 ;
  x1  = 3 ;
  y1  = 2 ;
  th1 = 0.314159265358979 ;
  k1  = -0.500000000000000 ;
  f1B = 0.154000000000000 ;
  */

  // test se basta un arco
  //curve.setup_G1( x0, y0, th0, x1, y1, th1 ) ;
  //valueType hatk0 = curve.getKappa();
  //valueType hatk1 = curve.getKappa()+curve.getSmax()*curve.getKappa_D() ;
  // check if hatk0 == k0   AND
  // check if hatk1 == k1
  //k0 = curve.getKappa() ;
  //k1 = curve.getKappa()+curve.getSmax()*curve.getKappa_D() ;

  // test con 2 archi  
  g2solve2arc.setup( x0, y0, th0, k0, x1, y1, th1, k1 ) ;
  bool converged = g2solve2arc.solve();
  std::cout << "\n\nTry Solution with TWO arcs\n";
  std::cout << "\n\nS0\n" << g2solve2arc.getS0();
  std::cout << "\n\nS1\n" << g2solve2arc.getS1();

  // test 3 archi
  g2solve3arc.setup( x0, y0, th0, k0, f0A, x1, y1, th1, k1, f1B ) ;
  converged = g2solve3arc.solve() ;
  
  Clothoid::ClothoidCurve const & S0 = g2solve3arc.getS0() ;
  Clothoid::ClothoidCurve const & S1 = g2solve3arc.getS1() ;
  Clothoid::ClothoidCurve const & SM = g2solve3arc.getSM() ;

  std::cout << "\n\nTry Solution with THREE arcs\n";
  std::cout << "\n\nS0 (NEW)\n" << S0 ;
  std::cout << "\n\nSM (NEW)\n" << SM ;
  std::cout << "\n\nS1 (NEW)\n" << S1 ;

  std::cout << "theta = " << S0.theta(S0.getSmax()) << " " << SM.theta(SM.getSmin()) << '\n' ;
  std::cout << "theta = " << SM.theta(SM.getSmax()) << " " << S1.theta(S1.getSmin()) << '\n' ;

  std::cout << "theta_D = " << S0.theta_D(S0.getSmax()) << " " << SM.theta_D(SM.getSmin()) << '\n' ;
  std::cout << "theta_D = " << SM.theta_D(SM.getSmax()) << " " << S1.theta_D(S1.getSmin()) << '\n' ;

  std::cout << "x = " << S0.X(S0.getSmax()) << " " << SM.X(SM.getSmin()) << '\n' ;
  std::cout << "y = " << SM.Y(SM.getSmax()) << " " << S1.Y(S1.getSmin()) << '\n' ;

  std::cout << "Dtheta = " << S0.theta(S0.getSmax())-SM.theta(SM.getSmin()) << '\n' ;
  std::cout << "Dtheta = " << SM.theta(SM.getSmax())-S1.theta(S1.getSmin()) << '\n' ;

  system("pause");

  return 0 ;
}
