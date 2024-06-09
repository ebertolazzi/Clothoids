//#define _USE_MATH_DEFINES
#include "Clothoids.hh"
#include "Clothoids_fmt.hh"

using G2lib::real_type;
using namespace std;

int
main() {

  G2lib::Biarc         ba{"temporary"};
  G2lib::ClothoidCurve c{"temporary"};

#if 0
  real_type x0  = 0;
  real_type y0  = 0;
  real_type th0 = Utils::m_pi/2;
  real_type x1  = 2;
  real_type y1  = 0;
  real_type th1 = -Utils::m_pi/2;
#else
  real_type x0  = -1;
  real_type y0  = 0;
  real_type th0 = Utils::m_pi/12;

  real_type x1  = 1;
  real_type y1  = 0;
  real_type th1 = -Utils::m_pi/4;
#endif

  ba.build( x0, y0, th0, x1, y1, th1 );
  c.build_G1( x0, y0, th0, x1, y1, th1 );

  cout << "BIARC\n";
  cout << ba;
  cout << "\n\n\nCLOTHOID\n";
  cout << c;

  cout << "\n\nALL DONE FOLKS!!!\n";

  return 0;
}
