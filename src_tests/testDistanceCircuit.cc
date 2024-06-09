//#define _USE_MATH_DEFINES
#include "Clothoids.hh"
#include "Clothoids_fmt.hh"

#include <sstream>

using G2lib::real_type;
using G2lib::integer;
using namespace std;

int
main() {

  try {

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

    //string fname = "circuit-fiorano_sim_kerbs_rebuilt.txt";
    string fname = "circuit-Fiorano_GoogleEarth_edges_no_kerbs_rebuilt.txt";
    ifstream file( fname.c_str() );

    UTILS_ASSERT( file.good(), "Cant find file: {}\n", fname );

    bool skipped_header = false;

    while ( file.good() ) {
      string str;
      getline(file,str);
      if ( str[0] == '#' ) {
        fmt::print( "SKIP LINE: {}\n", str );
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
      S.emplace_back(abscissa);
      K.emplace_back(curvature);
      X.emplace_back(x_mid_line);
      Y.emplace_back(y_mid_line);
    }

    G2lib::ClothoidList C{"temporary"};

    fmt::print( "x = {}\n", X[2230]-X[2231] );
    fmt::print( "y = {}\n", Y[2230]-Y[2231] );

    bool ok = C.build( x0, y0, theta0, S, K );

    fmt::print( "build = {}\n", ( ok ? "OK" : "NO OK" ) );

    #if 0
    real_type qx1 = 2.960934079000000e+02;
    real_type qy1 = 16.391055370000000;
    //real_type qx2 = 2.965952244000000e+02;
    //real_type qy2 = 16.475130419999999;

    integer idx = C.closest_segment( qx1, qy1 );
    fmt::print( "idx = {}\n", idx );

    real_type x, y, s, t, dst;
    integer  icurve;
    integer  icurve_begin = 140;
    integer  icurve_end   = 200;
    idx = C.closest_point_in_range_ISO(
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

    integer s_begin = 200;
    integer s_end   = 300;
    idx = C.closest_point_in_s_range_ISO(
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

    #endif

    // Nuovo test - - - - - - - - -
    real_type qx1 = 140.0; //100.0; //140.0;
    real_type qy1 = 10.0;  //-100; // 10.0;

    real_type x, y, s, t, dst;
    integer   icurve;
    integer   idx = C.closest_point_ISO(qx1, qy1, x, y, s, t, dst);

    fmt::print(
      "\n\n\nclosest_point_ISO for:({},{})\n"
      "idx    = {}\n"
      "x      = {}\n"
      "y      = {}\n"
      "s      = {}\n"
      "t      = {}\n"
      "dst    = {}\n",
      qx1,qy1,idx, x, y, s, t, dst, icurve
    );

    qx1 = 130.0;//80; //135;
    qy1 = -40;//-120.0;//30;
    real_type s_begin = s-10.0;
    real_type s_end   = s+100.0;
    idx = C.closest_point_in_s_range_ISO(
     qx1, qy1, s_begin, s_end, x, y, s, t, dst, icurve
    ) ;

    fmt::print(
      "closest_point_in_s_range_ISO for:({},{})\n"
      "idx    = {}\n"
      "x      = {}\n"
      "y      = {}\n"
      "s      = {}\n"
      "t      = {}\n"
      "dst    = {}\n"
      "icurve = {}\n",
      qx1, qy1, idx, x, y, s, t, dst, icurve
    );

  }
  catch ( std::exception const & e ) {
    std::cerr << e.what() << '\n';
  } catch (...) {
    std::cerr << "Unknown error!\n";
  }

  fmt::print( "\n\nALL DONE FOLKS!!!\n" );

  return 0;
}
