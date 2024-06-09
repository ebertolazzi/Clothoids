//#define _USE_MATH_DEFINES
#include "Clothoids.hh"
#include "Clothoids_fmt.hh"

using G2lib::real_type;
using G2lib::integer;
using namespace std;

int
main() {

  G2lib::ClothoidCurve C0{"C0"}, C1{"C1"};
  real_type x0     = 0;
  real_type y0     = 0;
  real_type theta0 = Utils::m_pi*0.7;
  real_type k0     = -0.1;
  real_type dk0    = 0.05;
  real_type L0     = 100;
  C0.build( x0, y0, theta0, k0, dk0, L0 );

  real_type x1     = -13;
  real_type y1     = 0.3;
  real_type theta1 = Utils::m_pi/4;
  real_type k1     = 0.1;
  real_type dk1    = -0.05;
  real_type L1     = 100;
  C1.build( x1, y1, theta1, k1, dk1, L1 );

  G2lib::LineSegment L{"L"};
  L.build( x0, y0, 2, 2 );

  G2lib::IntersectList ilist;
  C0.intersect( C1, ilist );
  C1.intersect( C0, ilist );
  G2lib::intersect( &L, &C1, ilist );
  G2lib::intersect( &C1, &L, ilist );

  ilist.clear();
  C0.intersect( C1, ilist );

  fmt::print( "L0 = {}\n", C0.length() );
  fmt::print( "L1 = {}\n", C1.length() );

  for ( size_t i = 0; i < ilist.size(); ++i )
    fmt::print( "s1[{}] = {}\n", i, ilist[i].first );

  for ( size_t i = 0; i < ilist.size(); ++i )
    fmt::print( "s2[{}] = {}\n", i, ilist[i].second );

  for ( size_t i = 0; i < ilist.size(); ++i )
    fmt::print( "x[{}] = {:10} y=[{}] = {:10}\n"
                "x[{}] = {:10} y=[{}] = {:10}\n\n",
                i, C0.X(ilist[i].first),
                i, C0.Y(ilist[i].first),
                i, C1.X(ilist[i].second),
                i, C1.Y(ilist[i].second) );

  fmt::print( "\n\nALL DONE FOLKS!!!\n" );

  return 0;
}
