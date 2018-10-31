//#define _USE_MATH_DEFINES
#include "Biarc.hh"
#include <cmath>
#include <iostream>

using G2lib::real_type;
using namespace std;

int
main() {

  G2lib::Biarc ba;

#if 0
  real_type x0  = 0;
  real_type y0  = 0;
  real_type th0 = G2lib::m_pi/2;
  real_type x1  = 2;
  real_type y1  = 0;
  real_type th1 = -G2lib::m_pi/2;
#else
  real_type x0  = -1;
  real_type y0  = 0;
  real_type th0 = G2lib::m_pi/12;

  real_type x1  = 1;
  real_type y1  = 0;
  real_type th1 = -G2lib::m_pi/4;
#endif

  ba.build( x0, y0, th0, x1, y1, th1 );

  cout << ba;

  cout << "\n\nALL DONE FOLKS!!!\n";

  return 0;
}
