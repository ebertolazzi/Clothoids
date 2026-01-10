// #define _USE_MATH_DEFINES
#include "Clothoids.hh"
#include "Clothoids_fmt.hh"
#include "Utils_eigen.hh"

using G2lib::integer;
using G2lib::real_type;
using namespace std;

int main()
{
  G2lib::ClothoidSplineG2 g2spline;

  G2lib::ClothoidList S{ "S" };
  ifstream            file( "G2_test.txt" );
  S.load( file );
  file.close();
  S.info( cout );

  integer   N{ 10 };
  real_type tol{ 0.5 };
  real_type X[]{ 0, 1, 2, 3, 4, 4 + tol, 5 + tol, 6 + tol, 7 + tol, 8 + tol };
  real_type Y[]{ 0, 0, 0, 0, 0, 1, 1, 1, 1, 1 };

  Utils::TicToc tictoc;

  Eigen::VectorXd THETA;
  THETA.resize( N );

  tictoc.tic();
  g2spline.build_PN( N, X, Y, THETA.data(), G2lib::TargetType::P8 );
  tictoc.toc();

  fmt::print( "build_G2 elapsed = {} [ms]\n", tictoc.elapsed_ms() );
  std::cout << THETA.transpose() << '\n';
  g2spline.info( std::cout );

  fmt::print( "\n\nALL DONE FOLKS!!!\n" );

  return 0;
}
