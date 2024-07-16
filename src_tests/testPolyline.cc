//#define _USE_MATH_DEFINES
#include "Clothoids.hh"
#include "Clothoids_fmt.hh"

using G2lib::real_type;
using namespace std;

int
main() {

  G2lib::ClothoidCurve C1{"C1"};
  G2lib::ClothoidCurve C2{"C2"};
  G2lib::PolyLine      P1{"P1"};
  G2lib::PolyLine      P2{"P2"};

  real_type x0  = 0;
  real_type y0  = 0;
  real_type th0 = 0;
  real_type x1  = 0;
  real_type y1  = 3;
  real_type th1 = Utils::m_pi;

  C1.build_G1( x0, y0, th0, x1, y1, th1 );
  C1.info(cout);

  C2 = C1;
  C2.rotate(Utils::m_pi/3,0,0);
  C2.translate(1,-1);

  G2lib::IntersectList ilist;
  C1.intersect( C2, ilist );

  fmt::print("CLOTHOIDS\n");

  for ( size_t i = 0; i < ilist.size(); ++i )
    fmt::print( "s1[{}] = {}\n", i, ilist[i].first );

  for ( size_t i = 0; i < ilist.size(); ++i )
    fmt::print( "s2[{}] = {}\n", i, ilist[i].second );

  for ( size_t i = 0; i < ilist.size(); ++i )
    fmt::print(
      "x[{}] = {} y[{}] = {}\n"
      "x[{}] = {} y[{}] = {}\n\n",
      i, C1.X(ilist[i].first),  i, C1.Y(ilist[i].first),
      i, C2.X(ilist[i].second), i, C2.Y(ilist[i].second)
    );

  P1.build( C1, 0.00001 );
  P2.build( C2, 0.00001 );

  fmt::print(
    "\n\nP2\n{}\n"
    "\n\nP1\n{}\n",
    P1.info(),
    P2.info()
  );

  vector<real_type> s1, s2;
  P1.intersect( P2, s1, s2 );

  fmt::print(
    "\n\nPOLYLINE\n"
    "# inter = {} {}\n",
    s1.size(), s2.size()
  );

  for ( auto is = s1.begin(); is != s1.end(); ++is )
    fmt::print( "s1[{}] = {}\n", is-s1.begin(), *is );

  for ( auto is = s2.begin(); is != s2.end(); ++is )
    fmt::print( "s2[{}] = {}\n", is-s2.begin(), *is );

  fmt::print("\n\nALL DONE FOLKS!!!\n");

  return 0;
}
