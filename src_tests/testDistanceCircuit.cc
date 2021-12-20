//#define _USE_MATH_DEFINES
#include "Clothoids.hh"
#include <sstream>

using G2lib::real_type;
using G2lib::int_type;
using namespace std;

int
main() {

  real_type x0     =   1.00793864;
  real_type y0     =   0.04992756;
  //real_type z0     =  -0.02035145;
  real_type theta0 =  -0.00034839;
  real_type x_mid_line, y_mid_line, dir_mid_line, abscissa,
            curvature, width_L, width_R, elevation, banking,
            slope, upsilon, torsion;

  vector<real_type> S;
  vector<real_type> K;
  vector<real_type> X;
  vector<real_type> Y;

  ifstream file("circuit-fiorano_sim_kerbs_rebuilt.txt");

  bool skipped_header = false;

  while ( file.good() ) {
    string str;
    getline(file,str);
    if ( str[0] == '#' ) {
      cout << "SKIP LINE: " << str << '\n';
      continue;
    }
    if ( !skipped_header ) { skipped_header = true; continue; }
    stringstream fstr(str);
    fstr
      >> x_mid_line
      >> y_mid_line
      >> dir_mid_line
      >> abscissa
      >> curvature
      >> width_L
      >> width_R
      >> elevation
      >> banking
      >> slope
      >> upsilon
      >> torsion;
    S.push_back(abscissa);
    K.push_back(curvature);
    X.push_back(x_mid_line);
    Y.push_back(y_mid_line);
  }

  G2lib::ClothoidList C;

  cout << "x = " << X[2230]-X[2231] << '\n';
  cout << "y = " << Y[2230]-Y[2231] << '\n';

  bool ok = C.build( x0, y0, theta0, S, K );

  cout << "build = " << ( ok ? "OK\n" : "NO OK\n" );

  real_type qx1 = 2.960934079000000e+02;
  real_type qy1 = 16.391055370000000;
  //real_type qx2 = 2.965952244000000e+02;
  //real_type qy2 = 16.475130419999999;

  int_type idx = C.closestSegment( qx1, qy1 );
  cout << "idx = " << idx << '\n';

  real_type x, y, s, t, dst;
  int_type  icurve;
  int_type  icurve_begin = 140;
  int_type  icurve_end   = 200;
  idx = C.closestPointInRange_ISO(
    qx1, qy1, icurve_begin, icurve_end, x, y, s, t, dst, icurve
  );

  fmt::print(
    "idx    = {}\n"
    "x      = {}\n"
    "y      = {}\n"
    "s      = {}\n"
    "t      = {}\n"
    "dst    = {}\n"
    "icurve = {}\n",
    idx, x, y, s, t, dst, icurve
  );

  int_type s_begin = 200;
  int_type s_end   = 300;
  idx = C.closestPointInSRange_ISO(
    qx1, qy1, s_begin, s_end, x, y, s, t, dst, icurve
  );

  fmt::print(
    "idx    = {}\n"
    "x      = {}\n"
    "y      = {}\n"
    "s      = {}\n"
    "t      = {}\n"
    "dst    = {}\n"
    "icurve = {}\n",
    idx, x, y, s, t, dst, icurve
  );

  cout << "\n\nALL DONE FOLKS!!!\n";

  return 0;
}
