//#define _USE_MATH_DEFINES
#include "Clothoids.hh"
#include "Clothoids_fmt.hh"
#include "Utils_eigen.hh"

using G2lib::real_type;
using G2lib::integer;
using namespace std;

int
main() {

  G2lib::G2solve2arc      g2solve2arc;
  G2lib::G2solve3arc      g2solve3arc;
  G2lib::ClothoidSplineG2 g2spline;

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
  
  integer N{27};
  real_type X[]{
     2.9265642,  2.6734362,  2.5109322,
     1.9078122,  1.1859282,  1.9249962,
     2.8265562,  0.0046842, -2.826567,
    -1.9437558, -1.1859438, -1.9062558,
    -2.501565,  -2.6734386, -2.9265642,
    -2.6187522, -1.1406318, -0.8968758,
    -1.4562558, -1.9062558, -0.0046878,
     1.9078122,  1.4468682,  0.8968722,
     1.1406282,  2.6187522,  2.9265642
  };
  real_type Y[]{
    -1.707808758, -1.707808758, -2.367185958,
    -2.582810358, -2.582810358, -1.167184758,
     0.915619242,  3.178123242,  0.915619242,
    -1.150000758, -2.582810358, -2.582810358,
    -2.393750358, -1.707808758, -1.707808758,
    -3.178123242, -3.178123242, -2.989063158,
    -0.915616758,  0.925003242,  2.953123242,
     0.925003242, -0.915616758, -2.989063158,
    -3.178123242, -3.178123242, -1.707808758
  };

  Eigen::VectorXd THETA;
  THETA.resize(27);

  g2spline.setP4();
  g2spline.build( N, X, Y, THETA.data() );
  
  std::cout << THETA.transpose() << '\n';
  
  g2spline.info( std::cout );

  fmt::print( "\n\nALL DONE FOLKS!!!\n" );

  return 0;
}
