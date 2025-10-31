//#define _USE_MATH_DEFINES
#include "Clothoids.hh"
#include "Clothoids_fmt.hh"
#include "Utils_eigen.hh"

using G2lib::real_type;
using G2lib::integer;
using namespace std;

int
main() {

  G2lib::ClothoidSplineG2 g2spline;

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
  
  Utils::TicToc tictoc;
  
  G2lib::TargetType targets[]{ G2lib::TargetType::P9,
                               G2lib::TargetType::P8,
                               G2lib::TargetType::P7,
                               G2lib::TargetType::P6,
                               G2lib::TargetType::P5,
                               G2lib::TargetType::P4 };

  //G2lib::TargetType targets[]{ G2lib::TargetType::P9 };

  for ( auto const & PN : targets ) {
  
    //auto PN{ G2lib::TargetType::P4 };

    Eigen::VectorXd THETA;
    THETA.resize(N);

    // Inizializza tutto a signaling NaN
    THETA.setConstant(std::numeric_limits<double>::signaling_NaN());

    tictoc.tic();
    g2spline.build_PN( N, X, Y, THETA.data(), PN );
    tictoc.toc();
    fmt::print( "build_G2 elapsed = {} [ms]\n", tictoc.elapsed_ms() );
    std::cout << THETA.transpose() << '\n';
    g2spline.info( std::cout );
  }

  fmt::print( "\n\nALL DONE FOLKS!!!\n" );

  return 0;
}

