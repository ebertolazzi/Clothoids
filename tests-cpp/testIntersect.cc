//#define _USE_MATH_DEFINES
#include "Clothoid.hh"
#include <cmath>
#include <iostream>

using G2lib::real_type;
using G2lib::int_type;

int
main() {

  G2lib::ClothoidCurve C0, C1;
  real_type x0     = -96.5947;
  real_type y0     = 88.6588;
  real_type theta0 = 1.65087;
  real_type k0     = 0.035746;
  real_type dk0    = 0;
  real_type L0     = 38.1914;
  C0.build( x0, y0, theta0, k0, dk0, L0 );

  real_type x1     = 0;
  real_type y1     = 0;
  real_type theta1 = 0;
  real_type k1     = 0.01;
  real_type dk1    = 0;
  real_type L1     = 600;
  C1.build( x1, y1, theta1, k1, dk1, L1 );

  std::vector<G2lib::real_type> s1, s2;
  G2lib::int_type  max_iter  = 100;
  G2lib::real_type tolerance = 1e-8;

  C0.intersect( C1, s1, s2, max_iter, tolerance );
  //C1.intersect( C0, s2, s1, max_iter, tolerance );

  std::cout << "L0 = " << C0.length() << '\n';
  std::cout << "L1 = " << C1.length() << '\n';

  for ( int_type i = 0; i < int_type(s1.size()); ++i )
    std::cout << "s1[ " << i << "] = " << s1[i] << '\n';

  for ( int_type i = 0; i < int_type(s2.size()); ++i )
    std::cout << "s2[ " << i << "] = " << s2[i] << '\n';

  for ( int_type i = 0; i < int_type(s1.size()); ++i )
    std::cout << "x["  << i << "] = " << C0.X(s1[i])
              << " y[" << i << "] = " << C0.Y(s1[i])
              << '\n';

  for ( int_type i = 0; i < int_type(s2.size()); ++i )
    std::cout << "x["  << i << "] = " << C1.X(s2[i])
              << " y[" << i << "] = " << C1.Y(s2[i])
              << '\n';

  std::cout << "\n\nALL DONE FOLKS!!!\n";

  return 0;
}
