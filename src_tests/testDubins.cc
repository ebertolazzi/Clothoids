//#define _USE_MATH_DEFINES
#include "Clothoids.hh"

using G2lib::real_type;
using G2lib::integer;
using namespace std;

int
main() {

  // check constructors
  real_type x0     = 0;
  real_type y0     = 0;
  real_type theta0 = 0*Utils::m_pi/2;
  real_type x3     = 1;
  real_type y3     = 0;
  real_type theta3 = -Utils::m_pi/2;
  real_type k_max  = 2;

  G2lib::Dubins DB;
  DB.build( x0, y0, theta0, x3, y3, theta3, k_max );
  fmt::print( "C0\n{}\n", DB.C0() );
  fmt::print( "C1\n{}\n", DB.C1() );
  fmt::print( "C2\n{}\n", DB.C2() );

  fmt::print( "\n\nALL DONE FOLKS!!!\n" );

  return 0;
}
