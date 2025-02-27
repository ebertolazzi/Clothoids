//#define _USE_MATH_DEFINES
#include "Clothoids.hh"
#include "Clothoids_fmt.hh"

using G2lib::real_type;
using G2lib::integer;
using namespace std;

int
main() {

  // check constructors
  constexpr real_type x0{0};
  constexpr real_type y0{2};
  constexpr real_type x1{4};
  constexpr real_type y1{3.5};
  constexpr real_type x2{6};
  constexpr real_type y2{1};

  G2lib::Triangle2D T1;
  T1.build( x0, y0, x1, y1, x2, y2, 0, 0, 0 );

  integer icode = T1.is_inside( T1.baricenter_x(), T1.baricenter_y() );

  fmt::print( "icode = {}\n", icode );

  icode = T1.is_inside( T1.P1() );

  fmt::print( "icode = {}\n", icode );

  icode = T1.is_inside( 10, 20 );

  fmt::print( "icode = {}\n", icode );

  {
    G2lib::ClothoidCurve C{"temporary"};

    constexpr real_type xx0    = 0;
    constexpr real_type yy0    = 2;
    constexpr real_type theta0 = 0;
    constexpr real_type kappa0 = 10;
    constexpr real_type dk     = -1;// 0;
    constexpr real_type L      = 1;
    real_type const max_angle { Utils::m_pi/18 };
    real_type const max_size  { 1e100 };
    C.build( xx0, yy0, theta0, kappa0, dk, L );

    vector<G2lib::Triangle2D> tvec;
    C.bb_triangles_ISO( 0, tvec, max_angle, max_size, 0 );

    for ( size_t i = 0; i < tvec.size(); ++i )
      fmt::print( "{} {}\n\n", i, tvec[i] );

    fmt::print( "{}\n",tvec.size() );
  }

  fmt::print( "\n\nALL DONE FOLKS!!!\n" );

  return 0;
}
