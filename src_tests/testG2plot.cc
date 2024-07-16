//#define _USE_MATH_DEFINES
#include "Clothoids.hh"
#include "Clothoids_fmt.hh"

#include <stack>
#include <ctime>

using G2lib::real_type;
using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846264338328
#endif

int
main() {

  G2lib::G2solve2arc g2sol;

  #if 1
  real_type x0  = -2;
  real_type y0  = 3;
  real_type th0 = M_PI/3;
  real_type k0  = 10.2;
  real_type x1  = 3;
  real_type y1  = 2;
  real_type th1 = M_PI/10;
  real_type k1  = -0.5;
  //real_type s1 = 1;
  #else
  real_type x0  = -1;
  real_type y0  = 0;
  real_type th0 = M_PI / 2;
  real_type k0  = -1;
  real_type s0  = 0.2;
  real_type x1  = 1;
  real_type y1  = 0;
  real_type th1 = -M_PI / 2;
  real_type k1  = -1;
  real_type s1  = 2;
  #endif

  //int iter = g2sol.build(x0, y0, th0, k0, s0, x1, y1, th1, k1, s1);
  int iter = g2sol.build(x0, y0, th0, k0, x1, y1, th1, k1 );
  //int iter = g2solve3arc.build( x0, y0, th0, k0, x1, y1, th1, k1 );
  fmt::print( "iter = {}\n", iter );

  G2lib::ClothoidCurve const & S0 = g2sol.S0();
  G2lib::ClothoidCurve const & S1 = g2sol.S1();

  fmt::print( "\n\nS0 (NEW)\n{}\n", S0 );
  fmt::print( "\n\nS1 (NEW)\n{}\n", S1 );

  fmt::print(
    "\n\n\n"
    "x1     = {}\n"
    "x0     = {}\n\n"
    "y1     = {}\n"
    "y0     = {}\n\n"
    "theta1 = {}\n"
    "theta0 = {}\n\n",
    S0.x_end(),
    S1.x_begin(),
    S0.y_end(),
    S1.y_begin(),
    S0.theta_end(),
    S1.theta_begin()
	);
	return 0;
}
