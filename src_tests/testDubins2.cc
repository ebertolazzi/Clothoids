//#define _USE_MATH_DEFINES
#include "Clothoids.hh"

using G2lib::real_type;
using G2lib::integer;
using namespace std;

#include <fstream>

int
main() {
  using std::abs;

  //ifstream file("DSdubinsP2P6Rand.csv");
  ifstream file("DSdubinsP2P6.csv");
  int wrong = 0;
  if (file.is_open()) {
    real_type x0, y0, th0, x1, y1, th1, kmax, l, t;
    std::string man;
    std::getline(file, man); // Discard the first line of the file

    real_type totTime = 0.0, totTime6man = 0.0;
    real_type minTime = std::numeric_limits<double>::max(), maxTime = 0.0;
    real_type minTime6man = std::numeric_limits<double>::max(), maxTime6man = 0.0;
    integer irow{0};
    while (file >> x0 >> y0 >> th0 >> x1 >> y1 >> th1 >> kmax >> l >> man >> t) {
      ++irow;
      bool error = false;

      auto startTime = std::chrono::high_resolution_clock::now();
      G2lib::Dubins dub(x0,y0, th0, x1, y1, th1, kmax);
      auto diff = std::chrono::duration<double, std::milli>(
            std::chrono::high_resolution_clock::now()-startTime
          ).count();
      totTime += diff;
      totTime6man += t;
      minTime = std::min(minTime, diff);
      maxTime = std::max(maxTime, diff);
      minTime6man = std::min(minTime6man, t);
      maxTime6man = std::max(maxTime6man, t);

      double totLen = dub.C0().length()+dub.C1().length()+dub.C2().length();

      char cman[3];
      cman[0] = (dub.C0().curvature()>0 ? 'L' : (dub.C0().curvature()<0 ? 'R' : 'S'));
      cman[1] = (dub.C1().curvature()>0 ? 'L' : (dub.C1().curvature()<0 ? 'R' : 'S'));
      cman[2] = (dub.C2().curvature()>0 ? 'L' : (dub.C2().curvature()<0 ? 'R' : 'S'));
      std::string compMan(cman, 3);

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
          "\nx0={} y0={} theta0={} x1={} y1={} theta1={} kmax={}\n",
          x0, y0, th0, x1, y1, th1, kmax
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
