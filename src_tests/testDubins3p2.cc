#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <string>
#include <iomanip>

#include "Clothoids.hh"
#include "Clothoids_fmt.hh"
#include "Utils_eigen.hh"

#include <thread>

using G2lib::real_type;
using G2lib::integer;
using namespace std;

int
main() {

  using std::abs;
  using Utils::m_pi;
  using Utils::m_2pi;

  // Sample one degree
  G2lib::Dubins3p DB_SD{"D3P_sample_one_degree"};
  DB_SD.set_sample_points(360);

  // Patern search
  G2lib::Dubins3p DB_PS{"D3P_pattern_search"};
  DB_PS.set_tolerance( 1e-10 );
  DB_PS.set_sample_angle( m_2pi/4 );

  G2lib::real_type x0     =                             -1;
  G2lib::real_type y0     =                              0;
  G2lib::real_type theta0 =         -1.4793802997313432179;
  G2lib::real_type xM     =        -0.94180274966005206316;
  G2lib::real_type yM     =        -0.94180274966005206316;
  G2lib::real_type xf     =                              1;
  G2lib::real_type yf     =                              0;
  G2lib::real_type thetaf =         -1.4793802997313432179;
  G2lib::real_type k_max  =         0.47036903761898174459;
  //G2lib::real_type thetaM =          4.5902159327450595683; (SAMPLE)
  //G2lib::real_type len    =          15.846598851367838634; (SAMPLE)
  //G2lib::real_type thetaM =          4.9447976995266156308; (PATTERN)
  //G2lib::real_type len    =          16.097284941403650294; (PATTERN)

  DB_SD.build( x0, y0, theta0, xM, yM, xf, yf, thetaf, k_max, G2lib::Dubins3pBuildType::SAMPLE_ONE_DEGREE );
  DB_PS.build( x0, y0, theta0, xM, yM, xf, yf, thetaf, k_max, G2lib::Dubins3pBuildType::PATTERN_SEARCH );

  std::cout << "SAMPLE\n";
  std::cout << DB_SD.info() << std::endl;
  std::cout << "PATTERN\n";
  std::cout << DB_PS.info() << std::endl;

  return 0;
}

