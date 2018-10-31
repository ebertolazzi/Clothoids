//#define _USE_MATH_DEFINES
#include "PolyLine.hh"
#include "Clothoid.hh"
#include <cmath>
#include <iostream>

using G2lib::real_type;
using namespace std;

int
main() {

  G2lib::ClothoidCurve C1;
  G2lib::ClothoidCurve C2;
  G2lib::PolyLine      P1;
  G2lib::PolyLine      P2;

  real_type x0  = 0;
  real_type y0  = 0;
  real_type th0 = 0;
  real_type x1  = 0;
  real_type y1  = 3;
  real_type th1 = G2lib::m_pi;

  C1.build_G1( x0, y0, th0, x1, y1, th1 );
  C1.info(cout);

  C2 = C1;
  C2.rotate(G2lib::m_pi/3,0,0);
  C2.translate(1,-1);

  vector<real_type> s1, s2;

  C1.intersect( C2, s1, s2, 30, 1e-10 );

  cout << "CLOTHOIDS\n";
  vector<real_type>::const_iterator is;
  for ( is = s1.begin(); is != s1.end(); ++is )
    cout << "s1[ " << is-s1.begin() << "] = "
         << *is << '\n';

  for ( is = s2.begin(); is != s2.end(); ++is )
    cout << "s2[ " << is-s2.begin() << "] = "
         << *is << '\n';

  P1.build( C1, 0.00001 );
  P2.build( C2, 0.00001 );

  cout << "\n\nP1\n";
  P1.info(cout);
  cout << "\n\nP2\n";
  P2.info(cout);

  s1.clear();
  s2.clear();
  P1.intersect( P2, s1, s2 );

  cout << "\n\nPOLYLINE\n";
  cout << "# inter = "
       << s1.size() << " "
       << s2.size() << '\n';

  for ( is = s1.begin(); is != s1.end(); ++is )
    cout << "s1[ " << is-s1.begin() << "] = "
         << *is << '\n';

  for ( is = s2.begin(); is != s2.end(); ++is )
    cout << "s2[ " << is-s2.begin() << "] = "
         << *is << '\n';

  G2lib::AABBtree aabb1, aabb2;
  P1.build_AABBtree( aabb1 );
  P2.build_AABBtree( aabb2 );

  //bool ok = aabb1.intersect(aabb2);
  //cout << "\nintersect = " << (ok?"YES\n":"NO\n");

  G2lib::AABBtree::VecPairPtrBBox intersectionList;
  aabb1.intersect( aabb2, intersectionList );
  cout << "# inter = " << intersectionList.size() << '\n';

  G2lib::AABBtree::VecPairPtrBBox::const_iterator ils;
  for ( ils = intersectionList.begin(); ils != intersectionList.end(); ++ils )
    cout
      << *(ils->first)
      << *(ils->second)
      << '\n';

  //P1.intersect(P2,s1,s2);
  //P2.info(cout);

  cout << "\n\nALL DONE FOLKS!!!\n";

  return 0;
}
