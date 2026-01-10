// #define _USE_MATH_DEFINES
#include "Clothoids.hh"
#include "Clothoids_fmt.hh"
#include "Utils_eigen.hh"

using G2lib::integer;
using G2lib::real_type;
using namespace std;

int main()
{
  G2lib::ClothoidList S{ "S" };

  integer   N{ 22 };
  real_type tol{ 0.5 };
  real_type x[]{ 0,        1,        2,        3,        4,        5,        6,        7,
                 8,        9,        10,       10 + tol, 11 + tol, 12 + tol, 13 + tol, 14 + tol,
                 15 + tol, 16 + tol, 17 + tol, 18 + tol, 19 + tol, 20 + tol };
  real_type y[]{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
  real_type w_min[]{ -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1,
                     -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1 };
  real_type w_max[]{ 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                     0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 };

  real_type theta0{ 0 };
  real_type theta1{ 0 };

  string target = "curvature";

  vector<real_type> theta( N );
  bool              ok = G2lib::build_guess_theta( N, x, y, theta.data() );

  Utils::TicToc tictoc;
  tictoc.tic();
  if ( ok )
  {
    if ( target == "lenght" )
    {
      auto target_fun = []( G2lib::ClothoidList const & lst ) -> real_type
      {
        real_type res{ 0 };
        for ( auto const & c : lst.get_list() ) res += c.length();
        return res;
      };
      ok = S.build_G2_with_target( N, x, y, theta.data(), w_min, w_max, theta0, theta1, target_fun );
    }
    else if ( target == "curvature" )
    {
      auto target_fun = []( G2lib::ClothoidList const & lst ) -> real_type
      {
        real_type res{ 0 };
        real_type len{ 0 };
        for ( auto const & c : lst.get_list() )
        {
          len += c.length();
          res += c.integral_curvature2();
        }
        return res / len;
      };
      ok = S.build_G2_with_target( N, x, y, theta.data(), w_min, w_max, theta0, theta1, target_fun );
    }
    else if ( target == "jerk" )
    {
      auto target_fun = []( G2lib::ClothoidList const & lst ) -> real_type
      {
        real_type res{ 0 };
        real_type len{ 0 };
        for ( auto const & c : lst.get_list() )
        {
          len += c.length();
          res += c.integral_jerk2();
        }
        return res / len;
      };
      ok = S.build_G2_with_target( N, x, y, theta.data(), w_min, w_max, theta0, theta1, target_fun );
    }
    else
    {
      UTILS_ASSERT0( false, "target must be in 'length', 'curvature', 'jerk'" );
    }
  }

  fmt::print( "build_G2 elapsed = {} [ms]\n", tictoc.elapsed_ms() );
  S.info( std::cout );

  fmt::print( "\n\nALL DONE FOLKS!!!\n" );

  return 0;
}
