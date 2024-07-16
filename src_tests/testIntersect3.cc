//#define _USE_MATH_DEFINES
#include "Clothoids.hh"
#include "Clothoids_fmt.hh"

using namespace std;
using G2lib::real_type;
using G2lib::integer;

int
main() {

  G2lib::LineSegment L0{"L0"};
  real_type x00     = 0;
  real_type y00     = 0;
  real_type theta00 = 3.14*0.9;
  L0.build( x00, y00, theta00, 10);

  G2lib::Biarc C0{"C0"}, C1{"C1"};
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

  G2lib::ClothoidList CL0(C0);
  G2lib::ClothoidList CL1(C1);

  G2lib::PolyLine PL0(C0,1e-8);
  G2lib::PolyLine PL1(C1,1e-4);

  G2lib::IntersectList ilist;

  G2lib::BaseCurve *pC0 = &PL0;
  G2lib::BaseCurve *pC1 = &PL1;
  //G2lib::BaseCurve *pC0 = &CL1;
  //G2lib::BaseCurve *pC1 = &L0;

  pC0->intersect( pC1,  ilist );
  pC0->intersect( &CL1, ilist );

  fmt::print( "{}\n", collision( pC0, pC1 ) );
  //C0.intersect( C1, ilist );

  fmt::print( "L0 = {}\n", C0.length() );
  fmt::print( "L1 = {}\n", C1.length() );

  for ( size_t i = 0; i < ilist.size(); ++i )
    fmt::print( "s1[ {} ] = {}\n", i, ilist[i].first );

  for ( size_t i = 0; i < ilist.size(); ++i )
    fmt::print( "s2[ {} ] = {}\n", i, ilist[i].second );

  for ( size_t i = 0; i < ilist.size(); ++i )
    fmt::print(
      "x[ {} ] = {:10}\n"
      "y[ {} ] = {:10}\n"
      "x[ {} ] = {:10}\n"
      "y[ {} ] = {:10}\n",
      i, C0.X(ilist[i].first),
      i, C0.Y(ilist[i].first),
      i, C1.X(ilist[i].second),
      i, C1.Y(ilist[i].second)
    );

  fmt::print( "\n\nALL DONE FOLKS!!!\n" );

  return 0;
}
