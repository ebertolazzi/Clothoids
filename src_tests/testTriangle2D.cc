//#define _USE_MATH_DEFINES
#include "Clothoids.hh"
#include "Clothoids_fmt.hh"

using G2lib::real_type;
using G2lib::integer;
using namespace std;

int
main() {

  // check constructors
  real_type x0{0};
  real_type y0{2};
  real_type x1{4};
  real_type y1{3.5};
  real_type x2{6};
  real_type y2{1};

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

    real_type xx0    = 0;
    real_type yy0    = 2;
    real_type theta0 = 0;
    real_type kappa0 = 10;
    real_type dk     = -1;// 0;
    real_type L      = 1;
    C.build( xx0, yy0, theta0, kappa0, dk, L );

    vector<G2lib::Triangle2D> tvec;
    C.bb_triangles_ISO( 0, tvec );

    for ( size_t i = 0; i < tvec.size(); ++i )
      fmt::print( "{} {}\n\n", i, tvec[i] );

    fmt::print( "{}\n",tvec.size() );
  }

  fmt::print( "\n\nALL DONE FOLKS!!!\n" );

  return 0;
}
