//#define _USE_MATH_DEFINES
#include "Clothoids.hh"
#include "Clothoids_fmt.hh"

using G2lib::real_type;
using G2lib::integer;

using namespace std;

int
main() {

  G2lib::ClothoidCurve C1{"temporary"};
  G2lib::ClothoidCurve C2{"temporary"};
  G2lib::ClothoidCurve C3{"temporary"};
  G2lib::ClothoidList  CL{"temporary"};
  {
  real_type x0     = 0.30002986753543353649;
  real_type y0     = -0.50753271613409067786;
  real_type theta0 = 1.7843235352254938064;
  real_type x1     = -2.9070958769989463377;
  real_type y1     = 14.283253009536736045;
  real_type theta1 = 1.8062988374013995152;
  C1.build_G1( x0, y0, theta0, x1, y1, theta1 );
  }
  {
  real_type x1     = -2.9070958769989481141;
  real_type y1     = 14.283253009536737821;
  real_type theta1 = 1.8062988374013995152;
  real_type x2     =  -5.1831205989265729528;
  real_type y2     = 22.926734844185681084;
  real_type theta2 = 1.8440574394333386632;
  C2.build_G1( x1, y1, theta1, x2, y2, theta2 );
  }
  {
  real_type x2     =  -5.183120598926573841;
  real_type y2     = 22.926734844185684636;
  real_type theta2 = 1.8440574394333386632;
  real_type x3     = -7.0388566648164294648;
  real_type y3     = 29.167179703253349743;
  real_type theta3 = 1.8598407392893718804;
  C3.build_G1( x2, y2, theta2, x3, y3, theta3 );
  }
  CL.push_back( C1 );
  CL.push_back( C2 );
  CL.push_back( C3 );

  std::cout << CL << "\n";

  for ( integer iii=0; iii<10; ++iii ) {
    CL.init();
    CL.push_back( C1 );
    CL.push_back( C2 );
    CL.push_back( C3 );

    real_type X, Y, S, T, DST;
    integer i = CL.closest_point_ISO( 0, 0, X, Y, S, T, DST );
    fmt::print(
      "i = {}\n"
      "X = {}\n"
      "Y = {}\n"
      "S = {}\n"
      "T = {}\n"
      "DST = {}\n",
      i, X, Y, S, T, DST
    );
  }

  fmt::print("\n\nALL DONE FOLKS!!!\n");

  return 0;
}
