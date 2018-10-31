//#define _USE_MATH_DEFINES
#include "Biarc.hh"
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;
using G2lib::real_type;
using G2lib::int_type;

int
main() {

  G2lib::Biarc C0, C1;
  real_type x0     = 0;
  real_type y0     = 0;
  real_type theta0 = 3.14*0.7;
  real_type x1     = 2;
  real_type y1     = 1;
  real_type theta1 = 3.14*0.1;
  C0.build( x0, y0, theta0, x1, y1, theta1 );

  x0     = -10;
  y0     = 0;
  theta0 = 3.14/4;
  x1     = 2;
  y1     = 1;
  theta1 = 3.14*0.1;
  C1.build( x0, y0, theta0, x1, y1, theta1 );

  G2lib::BaseCurve::IntersectList ilist;

  G2lib::BaseCurve *pC0 = &C0;
  G2lib::BaseCurve *pC1 = &C1;

  pC0->intersect( *pC1, ilist, false );
  //C0.intersect( C1, ilist );

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
