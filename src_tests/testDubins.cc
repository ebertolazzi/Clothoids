//#define _USE_MATH_DEFINES
#include "Clothoids.hh"
#include "Clothoids_fmt.hh"

using G2lib::real_type;
using G2lib::integer;
using namespace std;

int
main() {

  // check constructors
  constexpr real_type x0     = 484.405986676;
  constexpr real_type y0     = 556.309795113;
  constexpr real_type theta0 = 4.19036786001;
  constexpr real_type x3     = 486.116491569;
  constexpr real_type y3     = 556.286039501;
  constexpr real_type theta3 = 0.593139910671;
  constexpr real_type k_max  = 1;

  G2lib::Dubins DB{"DB"};
  bool ok { DB.build( x0, y0, theta0, x3, y3, theta3, k_max ) };
  fmt::print( "DB.build = {}\n", ok );
  fmt::print( "C0\n{}\n", DB.C0() );
  fmt::print( "C1\n{}\n", DB.C1() );
  fmt::print( "C2\n{}\n", DB.C2() );

  fmt::print( "\n\nALL DONE FOLKS!!!\n" );

  return 0;
}
