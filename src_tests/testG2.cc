//#define _USE_MATH_DEFINES
#include "Clothoids.hh"
#include "Clothoids_fmt.hh"

using G2lib::real_type;
using namespace std;

int
main() {

  G2lib::G2solve2arc g2solve2arc;
  G2lib::G2solve3arc g2solve3arc;

#if 1
  real_type x0  = -1;
  real_type y0  =  0;
  real_type th0 = -2;
  real_type k0  = 0.909297426825682;
  real_type x1  = 1;
  real_type y1  = 0;
  real_type th1 = 2;
  real_type k1  = 0.909297426825682;
#else
  real_type x0  = -1;
  real_type y0  = 0;
  real_type th0 = m_pi*0.9;
  real_type k0  = 0.2 + 1e-10;

  real_type x1  = 1;
  real_type y1  = 0;
  real_type th1 = -m_pi;
  real_type k1  = 0.2 + 0;
#endif

  // test 3 archi
  fmt::print( "\n\nTry Solution with THREE arcs\n" );
  int iter = g2solve3arc.build( x0, y0, th0, k0, x1, y1, th1, k1 );
  fmt::print( "iter = {}\n", iter );

  G2lib::ClothoidCurve const & S0 = g2solve3arc.S0();
  G2lib::ClothoidCurve const & S1 = g2solve3arc.S1();
  G2lib::ClothoidCurve const & SM = g2solve3arc.SM();

  fmt::print( "\n\nS0 (NEW)\n {}\n", S0 );
  fmt::print( "\n\nSM (NEW)\n {}\n", SM );
  fmt::print( "\n\nS1 (NEW)\n {}\n", S1 );

  fmt::print(
    "\n"
    "x  = {} {} err = {}\n"
    "y  = {} {} err = {}\n"
    "th = {} {} err = {}\n"
    "x  = {} {} err = {}\n"
    "y  = {} {} err = {}\n"
    "th = {} {} err = {}\n",
    S0.x_end(),       SM.x_begin(),     S0.x_end()-SM.x_begin(),
    S0.y_end(),       SM.y_begin(),     S0.y_end()-SM.y_begin(),
    S0.theta_end(),   SM.theta_begin(), S0.theta_end()-SM.theta_begin(),
    S1.x_begin(),     SM.x_end(),       S1.x_begin()-SM.x_end(),
    S1.y_begin(),     SM.y_end(),       S1.y_begin()-SM.y_end(),
    S1.theta_begin(), SM.theta_end(),   S1.theta_begin()-SM.theta_end()
  );

  G2lib::ClothoidList S{"S"};
  ifstream file("G2_test.txt");
  S.load(file);
  file.close();
  S.info( cout );

  fmt::print( "\n\nALL DONE FOLKS!!!\n" );

  return 0;
}
