//#define _USE_MATH_DEFINES
#include "Clothoids.hh"
#include "Clothoids_fmt.hh"

using G2lib::real_type;
using G2lib::integer;
using namespace std;

#include <fstream>

int
main() {
  using std::abs;
  using Utils::m_pi;

  G2lib::Dubins DB1{"DB1"};
  G2lib::Dubins DB2{"DB2"};

  {
    real_type k_max  = 1;
    real_type x0     = 0;
    real_type y0     = 0;
    real_type x3     = x0+1;
    real_type y3     = y0;
    real_type theta0 = -m_pi/2;
    real_type theta3 = m_pi/2;
    DB1.build( x0, y0, theta0, x3, y3, theta3, k_max );
  }
  {
    real_type k_max  = 1;
    real_type x0     = 1;
    real_type y0     = 0;
    real_type x3     = x0+1;
    real_type y3     = y0;
    real_type theta0 = m_pi/2;
    real_type theta3 = m_pi/2;
    DB2.build( x0, y0, theta0, x3, y3, theta3, k_max );
  }

  G2lib::IntersectList ilist;
  DB1.intersect( DB2, ilist );

  fmt::print( "{}\n", G2lib::collision( &DB1, &DB2 ) );

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
      i, DB1.X(ilist[i].first),
      i, DB1.Y(ilist[i].first),
      i, DB2.X(ilist[i].second),
      i, DB2.Y(ilist[i].second)
    );

  fmt::print( "\n\nALL DONE FOLKS!!!\n" );

  return 0;
}
