//#define _USE_MATH_DEFINES
#include "Clothoids.hh"

using G2lib::real_type;
using namespace std;

int
main() {

  G2lib::G2solve2arc g2solve2arc;
  G2lib::G2solve3arc g2solve3arc;

#if 1
  real_type x0  = -1;
  real_type y0  =  0;
  real_type th0 = -2;
  real_type k0  = 0.909297426825682;
  real_type x1  = 1;
  real_type y1  = 0;
  real_type th1 = 2;
  real_type k1  = 0.909297426825682;
#else
  real_type x0  = -1;
  real_type y0  = 0;
  real_type th0 = m_pi*0.9;
  real_type k0  = 0.2 + 1e-10;

  real_type x1  = 1;
  real_type y1  = 0;
  real_type th1 = -m_pi;
  real_type k1  = 0.2 + 0;
#endif

  // test 3 archi
  cout << "\n\nTry Solution with THREE arcs\n";
  int iter = g2solve3arc.build( x0, y0, th0, k0, x1, y1, th1, k1 );
  cout << "iter = " << iter << '\n';
  
  G2lib::ClothoidCurve const & S0 = g2solve3arc.getS0();
  G2lib::ClothoidCurve const & S1 = g2solve3arc.getS1();
  G2lib::ClothoidCurve const & SM = g2solve3arc.getSM();

  cout << "\n\nS0 (NEW)\n" << S0;
  cout << "\n\nSM (NEW)\n" << SM;
  cout << "\n\nS1 (NEW)\n" << S1;

  cout
    << "\nx  = " << S0.xEnd()       << " " << SM.xBegin()     << " err = " << S0.xEnd()-SM.xBegin()
    << "\ny  = " << S0.yEnd()       << " " << SM.yBegin()     << " err = " << S0.yEnd()-SM.yBegin()
    << "\nth = " << S0.thetaEnd()   << " " << SM.thetaBegin() << " err = " << S0.thetaEnd()-SM.thetaBegin()
    << "\nx  = " << S1.xBegin()     << " " << SM.xEnd()       << " err = " << S1.xBegin()-SM.xEnd()
    << "\ny  = " << S1.yBegin()     << " " << SM.yEnd()       << " err = " << S1.yBegin()-SM.yEnd()
    << "\nth = " << S1.thetaBegin() << " " << SM.thetaEnd()   << " err = " << S1.thetaBegin()-SM.thetaEnd()
    << '\n';

  G2lib::ClothoidList S;
  ifstream file("G2_test.txt");
  S.load(file);
  file.close();
  S.info( cout );

  cout << "\n\nALL DONE FOLKS!!!\n";

  return 0;
}
