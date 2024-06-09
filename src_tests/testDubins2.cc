//#define _USE_MATH_DEFINES
#include "Clothoids.hh"
#include "Clothoids_fmt.hh"

using G2lib::real_type;
using G2lib::integer;
using namespace std;

#include <fstream>

int
main() {
  using std::abs;

  ifstream file("DSdubinsP2P6Rand.csv");
  //ifstream file("DSdubinsP2P6.csv");
  int wrong = 0;
  if (file.is_open()) {
    real_type x0, y0, th0, x1, y1, th1, kmax, l, t;
    std::string man;
    std::getline(file, man); // Discard the first line of the file

    integer irow{0};
    while (file >> x0 >> y0 >> th0 >> x1 >> y1 >> th1 >> kmax >> l >> man >> t) {
      ++irow;
      bool error = false;

      G2lib::Dubins dub( x0, y0, th0, x1, y1, th1, kmax, "temporary" );
      double totLen = dub.length();
      std::string compMan(to_string(dub.solution_type()));

      if ( abs(totLen-l) > 1e-8) {
        fmt::print(
          "At row {} lengths are not the same: F/B = {}/{} diff = {}",
          irow, l, totLen, totLen-l
        );
        error = true;
      }
      if ( compMan != man ){
        fmt::print( " F:{}/B:{}", man, compMan );
        error = true;
      }
      if (error) {
        fmt::print(
          "\nx0={} y0={} theta0={} x1={} y1={} theta1={} kmax={} l={}:{}:{} errx={}, erry{}\n",
          x0, y0, th0, x1, y1, th1, kmax,
          dub.C0().length(), dub.C1().length(), dub.C2().length(),
          x1-dub.C2().x_end(), y1-dub.C2().y_end()
        );
        ++wrong;
      }

    }
    file.close();
    std::cout << "Wrong = " << wrong << std::endl;
  }
  else {
    std::cout << "Unable to open file" << std::endl;
  }
  return 0;
}
