//#define _USE_MATH_DEFINES
#include "Clothoid.hh"
#include <cmath>
#include <iostream>

using G2lib::real_type;
using G2lib::int_type;
using namespace std;

int
main() {

  G2lib::ClothoidCurve C0, C1;
  real_type x0     = 0;
  real_type y0     = 0;
  real_type theta0 = G2lib::m_pi*0.7;
  real_type k0     = -0.1;
  real_type dk0    = 0.05;
  real_type L0     = 100;
  C0.build( x0, y0, theta0, k0, dk0, L0 );

  real_type x1     = -13;
  real_type y1     = 0.3;
  real_type theta1 = G2lib::m_pi/4;
  real_type k1     = 0.1;
  real_type dk1    = -0.05;
  real_type L1     = 100;
  C1.build( x1, y1, theta1, k1, dk1, L1 );

  //vector<G2lib::real_type> s1, s2;
  //G2lib::int_type  max_iter  = 100;
  //G2lib::real_type tolerance = 1e-8;
  //C0.intersect( C1, s1, s2, max_iter, tolerance );
  //C1.intersect( C0, s2, s1, max_iter, tolerance );

  G2lib::IntersectList ilist;

  G2lib::noAABBtree();

  C0.intersect( C1, ilist, false );

  cout << "L0 = " << C0.length() << '\n';
  cout << "L1 = " << C1.length() << '\n';

  for ( size_t i = 0; i < ilist.size(); ++i )
    cout << "s1[ " << i << "] = " << ilist[i].first << '\n';

  for ( size_t i = 0; i < ilist.size(); ++i )
    cout << "s2[ " << i << "] = " << ilist[i].second << '\n';

  for ( size_t i = 0; i < ilist.size(); ++i )
    cout << "x[" << i << "] = " << setw(10) << C0.X(ilist[i].first)
         << " y[" << i << "] = " << setw(10) << C0.Y(ilist[i].first)
         << "\nx[" << i << "] = " << setw(10) << C1.X(ilist[i].second)
         << " y[" << i << "] = " << setw(10) << C1.Y(ilist[i].second)
         << "\n\n";

  cout << "\n\nALL DONE FOLKS!!!\n";

  return 0;
}
