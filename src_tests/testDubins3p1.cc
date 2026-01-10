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

using G2lib::integer;
using G2lib::real_type;
using namespace std;

#define PRINT_EVERY 1000

void printStats( int ith, Eigen::MatrixXd & times, Eigen::MatrixXd & ith_num, Eigen::MatrixXd & thetaM );
void printStatsFinal( Eigen::MatrixXd & times, Eigen::MatrixXd & ith_num, Eigen::MatrixXd & thetaM );

int main()
{
  using std::abs;
  using Utils::m_2pi;
  using Utils::m_pi;

  // Sample one degree
  G2lib::Dubins3p DB_SD{ "D3P_sample_one_degree" };
  DB_SD.set_sample_points( 360 );

  real_type tol{ 0.1 * m_pi / 180.0 };
  real_type sang{ m_2pi / 4 };

  // Patern search
  G2lib::Dubins3p DB_PS{ "D3P_pattern_search" };
  DB_PS.set_tolerance( tol );
  DB_PS.set_sample_angle( sang );
  G2lib::Dubins3p DB_PS_BRACKET{ "D3P_pattern_search_with_BRACKET" };
  DB_PS_BRACKET.set_tolerance( tol );
  DB_PS_BRACKET.set_sample_angle( sang );

  // Pattern trichotomy
  G2lib::Dubins3p DB_PT{ "D3P_pattern_trichotomy" };
  DB_PT.set_tolerance( tol );
  DB_PT.set_sample_angle( sang );
  G2lib::Dubins3p DB_PT_BRACKET{ "D3P_pattern_trichotomy_with_BRACKET" };
  DB_PT_BRACKET.set_tolerance( tol );
  DB_PT_BRACKET.set_sample_angle( sang );

  // Ellipse
  G2lib::Dubins3p DB_EL{ "D3P_ellipse" };

  Utils::TicToc tictoc;

  // Define constants and intervals
  G2lib::integer   numpts{ 10000 };  // Assuming a value for numpts
  G2lib::real_type xi{ -1 };
  G2lib::real_type yi{ 0 };
  G2lib::real_type xf{ 1 };
  G2lib::real_type yf{ 0 };

  G2lib::real_type k_max{ 1 };
  G2lib::real_type xm{ 0 };
  G2lib::real_type ym{ 0 };
  G2lib::real_type thi{ -m_pi / 2 };
  G2lib::real_type thf{ m_pi / 2 };

  G2lib::real_type elapsed_time1{ 0 };
  G2lib::real_type elapsed_time2{ 0 };
  G2lib::real_type elapsed_time3{ 0 };
  G2lib::real_type elapsed_time4{ 0 };
  G2lib::real_type elapsed_time5{ 0 };
  G2lib::real_type elapsed_time6{ 0 };

  Eigen::MatrixXd elapseTimeTable( numpts, 6 );

  Eigen::MatrixXd thetaMTable( numpts, 6 );
  Eigen::MatrixXd evaluationNumberTable( numpts, 6 );

  Eigen::MatrixXd type0Table( numpts, 6 );
  Eigen::MatrixXd type1Table( numpts, 6 );

  Eigen::MatrixXd lengthsTable_SD( numpts, 6 );
  Eigen::MatrixXd lengthsTable_PS( numpts, 6 );
  Eigen::MatrixXd lengthsTable_PT( numpts, 6 );
  Eigen::MatrixXd lengthsTable_PS_BRACKET( numpts, 6 );
  Eigen::MatrixXd lengthsTable_PT_BRACKET( numpts, 6 );
  Eigen::MatrixXd lengthsTable_EL( numpts, 6 );

  Eigen::Vector2d xm_interval( -2, 2 );
  Eigen::Vector2d ym_interval( -2, 2 );
  Eigen::Vector2d thi_interval( -m_pi, m_pi );
  Eigen::Vector2d thf_interval( -m_pi, m_pi );
  Eigen::Vector2d kmax_interval( 0.1, 1.5 );

  // Seed the random number generator
  std::mt19937                              gen( 1 );
  std::uniform_real_distribution<real_type> dis( 0.0, 1.0 );

  // Create a matrix to store the data
  Eigen::MatrixXd dataTable( numpts, 9 );

  // Fill the matrix with data
  for ( int i{ 0 }; i < numpts; ++i )
  {
    G2lib::real_type ss{ static_cast<G2lib::real_type>( i ) / static_cast<G2lib::real_type>( numpts ) };
    // if (i % PRINT_EVERY == 0)  std::cout << "completed " << ( (double) i ) / ((double) numpts) << " %, " << "i = " <<
    // i << std::endl;
    if ( i % PRINT_EVERY == 0 ) fmt::print( "completed {:>5.3f} %, i = {}\n", 100 * ss, i );

    dataTable( i, 0 ) = xi;  // Pi_x
    dataTable( i, 1 ) = yi;  // Pi_y
    dataTable( i, 4 ) = xf;  // Pf_x
    dataTable( i, 5 ) = yf;  // Pf_y

    G2lib::real_type r{ dis( gen ) };

    thi   = thi_interval[0] + ( thi_interval[1] - thi_interval[0] ) * r;     // Pi_th
    xm    = xm_interval[0] + ( xm_interval[1] - xm_interval[0] ) * r;        // Pm_x
    ym    = ym_interval[0] + ( ym_interval[1] - ym_interval[0] ) * r;        // Pm_y
    thf   = thf_interval[0] + ( thf_interval[1] - thf_interval[0] ) * r;     // Pf_th
    k_max = kmax_interval[0] + ( kmax_interval[1] - kmax_interval[0] ) * r;  // kmax

    dataTable( i, 2 ) = thi;    // Pi_th
    dataTable( i, 3 ) = xm;     // Pm_x
    dataTable( i, 4 ) = ym;     // Pm_y
    dataTable( i, 6 ) = thf;    // Pf_th
    dataTable( i, 7 ) = k_max;  // kmax

    // ███████╗ █████╗ ███╗   ███╗██████╗ ██╗     ███████╗
    // ██╔════╝██╔══██╗████╗ ████║██╔══██╗██║     ██╔════╝
    // ███████╗███████║██╔████╔██║██████╔╝██║     █████╗
    // ╚════██║██╔══██║██║╚██╔╝██║██╔═══╝ ██║     ██╔══╝
    // ███████║██║  ██║██║ ╚═╝ ██║██║     ███████╗███████╗
    // ╚══════╝╚═╝  ╚═╝╚═╝     ╚═╝╚═╝     ╚══════╝╚══════╝

    // threads.emplace_back([&](){
    tictoc.tic();
    DB_SD.build( xi, yi, thi, xm, ym, xf, yf, thf, k_max, G2lib::Dubins3pBuildType::SAMPLE_ONE_DEGREE );
    tictoc.toc();
    elapsed_time1 = tictoc.elapsed_ms();

    lengthsTable_SD( i, 0 ) = DB_SD.length0();
    lengthsTable_SD( i, 1 ) = DB_SD.length1();
    lengthsTable_SD( i, 2 ) = DB_SD.length2();
    lengthsTable_SD( i, 3 ) = DB_SD.length3();
    lengthsTable_SD( i, 4 ) = DB_SD.length4();
    lengthsTable_SD( i, 5 ) = DB_SD.length5();

    // ██████╗  █████╗ ████████╗████████╗███████╗██████╗ ███╗   ██╗    ███████╗███████╗ █████╗ ██████╗  ██████╗██╗  ██╗
    // ██╔══██╗██╔══██╗╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗████╗  ██║    ██╔════╝██╔════╝██╔══██╗██╔══██╗██╔════╝██║  ██║
    // ██████╔╝███████║   ██║      ██║   █████╗  ██████╔╝██╔██╗ ██║    ███████╗█████╗  ███████║██████╔╝██║     ███████║
    // ██╔═══╝ ██╔══██║   ██║      ██║   ██╔══╝  ██╔══██╗██║╚██╗██║    ╚════██║██╔══╝  ██╔══██║██╔══██╗██║     ██╔══██║
    // ██║     ██║  ██║   ██║      ██║   ███████╗██║  ██║██║ ╚████║    ███████║███████╗██║  ██║██║  ██║╚██████╗██║  ██║
    // ╚═╝     ╚═╝  ╚═╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚═╝  ╚═══╝    ╚══════╝╚══════╝╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝

    tictoc.tic();
    DB_PS.build( xi, yi, thi, xm, ym, xf, yf, thf, k_max, G2lib::Dubins3pBuildType::PATTERN_SEARCH );
    tictoc.toc();
    elapsed_time2 = tictoc.elapsed_ms();

    lengthsTable_PS( i, 0 ) = DB_PS.length0();
    lengthsTable_PS( i, 1 ) = DB_PS.length1();
    lengthsTable_PS( i, 2 ) = DB_PS.length2();
    lengthsTable_PS( i, 3 ) = DB_PS.length3();
    lengthsTable_PS( i, 4 ) = DB_PS.length4();
    lengthsTable_PS( i, 5 ) = DB_PS.length5();

    tictoc.tic();
    DB_PS_BRACKET
      .build( xi, yi, thi, xm, ym, xf, yf, thf, k_max, G2lib::Dubins3pBuildType::PATTERN_SEARCH_WITH_ALGO_BRACKET );
    tictoc.toc();
    elapsed_time3 = tictoc.elapsed_ms();

    lengthsTable_PS_BRACKET( i, 0 ) = DB_PS_BRACKET.length0();
    lengthsTable_PS_BRACKET( i, 1 ) = DB_PS_BRACKET.length1();
    lengthsTable_PS_BRACKET( i, 2 ) = DB_PS_BRACKET.length2();
    lengthsTable_PS_BRACKET( i, 3 ) = DB_PS_BRACKET.length3();
    lengthsTable_PS_BRACKET( i, 4 ) = DB_PS_BRACKET.length4();
    lengthsTable_PS_BRACKET( i, 5 ) = DB_PS_BRACKET.length5();

    // ██████╗  █████╗ ████████╗████████╗███████╗██████╗ ███╗   ██╗    ████████╗██████╗ ██╗ ██████╗
    // ██╔══██╗██╔══██╗╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗████╗  ██║    ╚══██╔══╝██╔══██╗██║██╔════╝
    // ██████╔╝███████║   ██║      ██║   █████╗  ██████╔╝██╔██╗ ██║       ██║   ██████╔╝██║██║
    // ██╔═══╝ ██╔══██║   ██║      ██║   ██╔══╝  ██╔══██╗██║╚██╗██║       ██║   ██╔══██╗██║██║
    // ██║     ██║  ██║   ██║      ██║   ███████╗██║  ██║██║ ╚████║       ██║   ██║  ██║██║╚██████╗
    // ╚═╝     ╚═╝  ╚═╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚═╝  ╚═══╝       ╚═╝   ╚═╝  ╚═╝╚═╝ ╚═════╝

    tictoc.tic();
    DB_PT.build( xi, yi, thi, xm, ym, xf, yf, thf, k_max, G2lib::Dubins3pBuildType::PATTERN_TRICHOTOMY );
    tictoc.toc();
    elapsed_time4 = tictoc.elapsed_ms();

    lengthsTable_PT( i, 0 ) = DB_PT.length0();
    lengthsTable_PT( i, 1 ) = DB_PT.length1();
    lengthsTable_PT( i, 2 ) = DB_PT.length2();
    lengthsTable_PT( i, 3 ) = DB_PT.length3();
    lengthsTable_PT( i, 4 ) = DB_PT.length4();
    lengthsTable_PT( i, 5 ) = DB_PT.length5();

    tictoc.tic();
    DB_PT_BRACKET
      .build( xi, yi, thi, xm, ym, xf, yf, thf, k_max, G2lib::Dubins3pBuildType::PATTERN_TRICHOTOMY_WITH_ALGO_BRACKET );
    tictoc.toc();
    elapsed_time5 = tictoc.elapsed_ms();

    lengthsTable_PT_BRACKET( i, 0 ) = DB_PT_BRACKET.length0();
    lengthsTable_PT_BRACKET( i, 1 ) = DB_PT_BRACKET.length1();
    lengthsTable_PT_BRACKET( i, 2 ) = DB_PT_BRACKET.length2();
    lengthsTable_PT_BRACKET( i, 3 ) = DB_PT_BRACKET.length3();
    lengthsTable_PT_BRACKET( i, 4 ) = DB_PT_BRACKET.length4();
    lengthsTable_PT_BRACKET( i, 5 ) = DB_PT_BRACKET.length5();

    if ( 1.0014 * DB_SD.length() < DB_PT.length() )
    {
      fmt::print(
        "x₀   = {:30.20};\n"
        "y₀   = {:30.20};\n"
        "θ₀   = {:30.20};\n"
        "xᵐ   = {:30.20};\n"
        "yᵐ   = {:30.20};\n"
        "xᶠ   = {:30.20};\n"
        "yᶠ   = {:30.20};\n"
        "θf   = {:30.20};\n"
        "κᵐᵃˣ = {:30.20};\n"
        "θᵐ   = {:30.20}; (SAMPLE)\n"
        "len  = {:30.20}; (SAMPLE)\n"
        "θᵐ   = {:30.20}; (PATTERN)\n"
        "len  = {:30.20}; (PATTERN)\n\n",
        xi,
        yi,
        thi,
        xm,
        ym,
        xf,
        yf,
        thf,
        k_max,
        DB_SD.theta3_begin(),
        DB_SD.length(),
        DB_PT.theta3_begin(),
        DB_PT.length() );
    }

    // ███████╗██╗     ██╗     ██╗██████╗ ███████╗███████╗
    // ██╔════╝██║     ██║     ██║██╔══██╗██╔════╝██╔════╝
    // █████╗  ██║     ██║     ██║██████╔╝███████╗█████╗
    // ██╔══╝  ██║     ██║     ██║██╔═══╝ ╚════██║██╔══╝
    // ███████╗███████╗███████╗██║██║     ███████║███████╗
    // ╚══════╝╚══════╝╚══════╝╚═╝╚═╝     ╚══════╝╚══════╝

    tictoc.tic();
    DB_EL.build( xi, yi, thi, xm, ym, xf, yf, thf, k_max, G2lib::Dubins3pBuildType::ELLIPSE );
    tictoc.toc();
    elapsed_time6 = tictoc.elapsed_ms();

    lengthsTable_EL( i, 0 ) = DB_EL.length0();
    lengthsTable_EL( i, 1 ) = DB_EL.length1();
    lengthsTable_EL( i, 2 ) = DB_EL.length2();
    lengthsTable_EL( i, 3 ) = DB_EL.length3();
    lengthsTable_EL( i, 4 ) = DB_EL.length4();
    lengthsTable_EL( i, 5 ) = DB_EL.length5();

    // ███████╗███╗   ██╗███████╗███████╗███╗   ███╗██████╗ ██╗     ███████╗
    // ██╔════╝████╗  ██║██╔════╝██╔════╝████╗ ████║██╔══██╗██║     ██╔════╝
    // █████╗  ██╔██╗ ██║███████╗█████╗  ██╔████╔██║██████╔╝██║     █████╗
    // ██╔══╝  ██║╚██╗██║╚════██║██╔══╝  ██║╚██╔╝██║██╔══██╗██║     ██╔══╝
    // ███████╗██║ ╚████║███████║███████╗██║ ╚═╝ ██║██████╔╝███████╗███████╗
    // ╚══════╝╚═╝  ╚═══╝╚══════╝╚══════╝╚═╝     ╚═╝╚═════╝ ╚══════╝╚══════╝

    elapseTimeTable( i, 0 ) = elapsed_time1;
    elapseTimeTable( i, 1 ) = elapsed_time2;
    elapseTimeTable( i, 2 ) = elapsed_time3;
    elapseTimeTable( i, 3 ) = elapsed_time4;
    elapseTimeTable( i, 4 ) = elapsed_time5;
    elapseTimeTable( i, 5 ) = elapsed_time6;

    thetaMTable( i, 0 ) = DB_SD.theta2( 0 );
    thetaMTable( i, 1 ) = DB_PS.theta2( 0 );
    thetaMTable( i, 2 ) = DB_PS_BRACKET.theta2( 0 );
    thetaMTable( i, 3 ) = DB_PT.theta2( 0 );
    thetaMTable( i, 4 ) = DB_PT_BRACKET.theta2( 0 );
    thetaMTable( i, 5 ) = DB_EL.theta2( 0 );

    evaluationNumberTable( i, 0 ) = DB_SD.num_evaluation();
    evaluationNumberTable( i, 1 ) = DB_PS.num_evaluation();
    evaluationNumberTable( i, 2 ) = DB_PS_BRACKET.num_evaluation();
    evaluationNumberTable( i, 3 ) = DB_PT.num_evaluation();
    evaluationNumberTable( i, 4 ) = DB_PT_BRACKET.num_evaluation();
    evaluationNumberTable( i, 5 ) = DB_EL.num_evaluation();

    type0Table( i, 0 ) = static_cast<double>( DB_SD.solution_type0() );
    type0Table( i, 1 ) = static_cast<double>( DB_PS.solution_type0() );
    type0Table( i, 2 ) = static_cast<double>( DB_PS_BRACKET.solution_type0() );
    type0Table( i, 3 ) = static_cast<double>( DB_PT.solution_type0() );
    type0Table( i, 4 ) = static_cast<double>( DB_PT_BRACKET.solution_type0() );
    type0Table( i, 5 ) = static_cast<double>( DB_EL.solution_type0() );

    type1Table( i, 0 ) = static_cast<double>( DB_SD.solution_type1() );
    type1Table( i, 1 ) = static_cast<double>( DB_PS.solution_type1() );
    type1Table( i, 2 ) = static_cast<double>( DB_PS_BRACKET.solution_type1() );
    type1Table( i, 3 ) = static_cast<double>( DB_PT.solution_type1() );
    type1Table( i, 4 ) = static_cast<double>( DB_PT_BRACKET.solution_type1() );
    type1Table( i, 5 ) = static_cast<double>( DB_EL.solution_type1() );

    printStats( i, elapseTimeTable, evaluationNumberTable, thetaMTable );
  }

  // Statistics
  printStatsFinal( elapseTimeTable, evaluationNumberTable, thetaMTable );

  Eigen::MatrixXd L_SD( numpts, 1 );
  Eigen::MatrixXd L_PS( numpts, 1 );
  Eigen::MatrixXd L_PS_BRACKET( numpts, 1 );
  Eigen::MatrixXd L_PT( numpts, 1 );
  Eigen::MatrixXd L_PT_BRACKET( numpts, 1 );
  Eigen::MatrixXd L_EL( numpts, 1 );
  L_SD         = lengthsTable_SD.rowwise().sum();
  L_PS         = lengthsTable_PS.rowwise().sum();
  L_PS_BRACKET = lengthsTable_PS_BRACKET.rowwise().sum();
  L_PT         = lengthsTable_PT.rowwise().sum();
  L_PT_BRACKET = lengthsTable_PT_BRACKET.rowwise().sum();
  L_EL         = lengthsTable_EL.rowwise().sum();

  Eigen::MatrixXd rappEL( numpts, 1 );
  Eigen::MatrixXd rappPS( numpts, 1 );
  Eigen::MatrixXd rappPS_BRACKET( numpts, 1 );
  Eigen::MatrixXd rappPT( numpts, 1 );
  Eigen::MatrixXd rappPT_BRACKET( numpts, 1 );

  rappEL         = L_EL.array() / L_SD.array();
  rappPS         = L_PS.array() / L_SD.array();
  rappPS_BRACKET = L_PS_BRACKET.array() / L_SD.array();
  rappPT         = L_PT.array() / L_SD.array();
  rappPT_BRACKET = L_PT_BRACKET.array() / L_SD.array();

  fmt::print(
    "L_SD = {:>10.5f}, "
    "L_PS = {:>10.5f}, "
    "L_PS_BRACKET = {:>10.5f}, "
    "L_PT = {:>10.5f}, "
    "L_PT_BRACKET = {:>10.5f}, "
    "L_EL = {:>10.5f}\n",
    L_SD.mean(),
    L_PS.mean(),
    L_PS_BRACKET.mean(),
    L_PT.mean(),
    L_PT_BRACKET.mean(),
    L_EL.mean() );
  // print max and min
  fmt::print(
    "L_SD = {:>10.5f}, "
    "L_PS = {:>10.5f}, "
    "L_PS_BRACKET = {:>10.5f}, "
    "L_PT = {:>10.5f}, "
    "L_PT_BRACKET = {:>10.5f}, "
    "L_EL = {:>10.5f}\n",
    L_SD.maxCoeff(),
    L_PS.maxCoeff(),
    L_PS_BRACKET.maxCoeff(),
    L_PT.maxCoeff(),
    L_PT_BRACKET.maxCoeff(),
    L_EL.maxCoeff() );
  fmt::print(
    "L_SD = {:>10.5f}, "
    "L_PS = {:>10.5f}, "
    "L_PS_BRACKET = {:>10.5f}, "
    "L_PT = {:>10.5f}, "
    "L_PT_BRACKET = {:>10.5f}, "
    "L_EL = {:>10.5f}\n",
    L_SD.minCoeff(),
    L_PS.minCoeff(),
    L_PS_BRACKET.minCoeff(),
    L_PT.minCoeff(),
    L_PT_BRACKET.minCoeff(),
    L_EL.minCoeff() );

  fmt::print(
    "rapp EL = {:>10.5f}:mean, {:>10.5f}:max, {:>10.5f}:min \n",
    rappEL.mean(),
    rappEL.maxCoeff(),
    rappEL.minCoeff() );

  fmt::print(
    "rapp PS = {:>10.5f}:mean, {:>10.5f}:max, {:>10.5f}:min \n",
    rappPS.mean(),
    rappPS.maxCoeff(),
    rappPS.minCoeff() );
  fmt::print(
    "rapp PS_BRACKET = {:>10.5f}:mean, {:>10.5f}:max, {:>10.5f}:min \n",
    rappPS_BRACKET.mean(),
    rappPS_BRACKET.maxCoeff(),
    rappPS_BRACKET.minCoeff() );

  fmt::print(
    "rapp PT = {:>10.5f}:mean, {:>10.5f}:max, {:>10.5f}:min \n",
    rappPT.mean(),
    rappPT.maxCoeff(),
    rappPT.minCoeff() );

  fmt::print(
    "rapp PT_BRACKET = {:>10.5f}:mean, {:>10.5f}:max, {:>10.5f}:min \n",
    rappPT_BRACKET.mean(),
    rappPT_BRACKET.maxCoeff(),
    rappPT_BRACKET.minCoeff() );

  G2lib::real_type dr{ 1e-4 };
  for ( real_type r{ 1.0 }; r < 10; r += dr )
  {
    integer o_EL{ static_cast<integer>( ( rappEL.array() > r ).count() ) };
    integer o_PS{ static_cast<integer>( ( rappPS.array() > r ).count() ) };
    integer o_PS_BRACKET{ static_cast<integer>( ( rappPS_BRACKET.array() > r ).count() ) };
    integer o_PT{ static_cast<integer>( ( rappPT.array() > r ).count() ) };
    integer o_PT_BRACKET{ static_cast<integer>( ( rappPT_BRACKET.array() > r ).count() ) };
    fmt::print(
      "num outliers > {:>12.8f}  "
      "EL = {:>6}  "
      "PS = {:>6}  "
      "PS_BRACKET = {:>6}  "
      "PT = {:>6}  "
      "PT_BRACKET = {:>6}\n",
      r,
      o_EL,
      o_PS,
      o_PS_BRACKET,
      o_PT,
      o_PT_BRACKET );
    if ( o_EL == 0 && o_PS == 0 && o_PS_BRACKET == 0 && o_PT == 0 && o_PT_BRACKET == 0 ) break;
    dr *= 2;
  }

  // Save the matrix to a CSV file
  std::ofstream file( "store_possible_combination.csv" );
  file << "Pi_x,Pi_y,Pi_th,Pm_x,Pm_y,Pf_x,Pf_y,Pf_th,kmax\n";
  for ( int i{ 0 }; i < numpts; ++i )
  {
    for ( int j{ 0 }; j < 9; ++j )
    {
      file << dataTable( i, j );
      if ( j < 8 ) file << ",";
    }
    file << "\n";
  }

  file.close();
  std::cout << "Data saved to store_possible_combination.csv" << std::endl;

  std::ofstream file_SD( "store_solution_SD.csv" );
  file_SD << "length0,length1,length2,length3,length4,length5,theta2,eval,type0,type1,time\n";
  for ( int i{ 0 }; i < numpts; ++i )
  {
    for ( int j{ 0 }; j < 6; ++j )
    {
      file_SD << lengthsTable_SD( i, j );
      if ( j < 5 ) file_SD << ",";
    }
    file_SD << "," << thetaMTable( i, 0 ) << "," << evaluationNumberTable( i, 0 ) << "," << type0Table( i, 0 ) << ","
            << type1Table( i, 0 ) << "," << elapseTimeTable( i, 0 ) << "\n";
  }

  file_SD.close();

  std::cout << "Data saved to store_solution_SD.csv" << std::endl;

  std::ofstream file_PS( "store_solution_PS.csv" );
  file_PS << "length0,length1,length2,length3,length4,length5,theta2,eval,type0,type1,time\n";
  for ( int i{ 0 }; i < numpts; ++i )
  {
    for ( int j{ 0 }; j < 6; ++j )
    {
      file_PS << lengthsTable_PS( i, j );
      if ( j < 5 ) file_PS << ",";
    }
    file_PS << "," << thetaMTable( i, 1 ) << "," << evaluationNumberTable( i, 1 ) << "," << type0Table( i, 1 ) << ","
            << type1Table( i, 1 ) << "," << elapseTimeTable( i, 1 ) << "\n";
  }
  file_PS.close();

  std::cout << "Data saved to store_solution_PS.csv" << std::endl;

  std::ofstream file_PT( "store_solution_PT.csv" );
  file_PT << "length0,length1,length2,length3,length4,length5,theta2,eval,type0,type1,time\n";
  for ( int i{ 0 }; i < numpts; ++i )
  {
    for ( int j{ 0 }; j < 6; ++j )
    {
      file_PT << lengthsTable_PT( i, j );
      if ( j < 5 ) file_PT << ",";
    }
    file_PT << "," << thetaMTable( i, 2 ) << "," << evaluationNumberTable( i, 2 ) << "," << type0Table( i, 2 ) << ","
            << type1Table( i, 2 ) << "," << elapseTimeTable( i, 2 ) << "\n";
  }

  file_PT.close();

  std::cout << "Data saved to store_solution_PT.csv" << std::endl;

  std::ofstream file_EL( "store_solution_EL.csv" );
  file_EL << "length0,length1,length2,length3,length4,length5,theta2,eval,type0,type1,time\n";
  for ( int i{ 0 }; i < numpts; ++i )
  {
    for ( int j{ 0 }; j < 6; ++j )
    {
      file_EL << lengthsTable_EL( i, j );
      if ( j < 5 ) file_EL << ",";
    }
    file_EL << "," << thetaMTable( i, 3 ) << "," << evaluationNumberTable( i, 3 ) << "," << type0Table( i, 3 ) << ","
            << type1Table( i, 3 ) << "," << elapseTimeTable( i, 3 ) << "\n";
  }
  file_EL.close();

  std::cout << "Data saved to store_solution_EL.csv" << std::endl;

  return 0;
}

void printStats( int ith, Eigen::MatrixXd & times, Eigen::MatrixXd & ith_num, Eigen::MatrixXd & thetaM )
{
  if ( ith % PRINT_EVERY != 0 ) return;
  fmt::print( "---------------------------------------------\n" );
  fmt::print( "Result stats iter: {:>10}:\n", ith );
  fmt::print(
    "  {:<10} {:<20} {:<20} {:<20} {:<20} {:<20} {:<20}\n",
    "VAL",
    "SD ",
    "PS ",
    "PS_BRACKET ",
    "PT ",
    "PT_BRACKET ",
    "ELL" );
  fmt::print( "  {:<10} {:<20} {:<20} {:<20} {:<20} {:<20} {:<20}\n", "---", "---", "---", "---", "---", "---", "---" );
  fmt::print(
    "  {:<10} {:<20.15g} {:<20.15g} {:<20.15g} {:<20.15g} {:<20.15g} {:<20.15g}\n",
    "time",
    times( ith, 0 ),
    times( ith, 1 ),
    times( ith, 2 ),
    times( ith, 3 ),
    times( ith, 4 ),
    times( ith, 5 ) );
  fmt::print(
    "  {:<10} {:<20.15g} {:<20.15g} {:<20.15g} {:<20.15g} {:<20.15g} {:<20.15g}\n",
    "eval",
    ith_num( ith, 0 ),
    ith_num( ith, 1 ),
    ith_num( ith, 2 ),
    ith_num( ith, 3 ),
    ith_num( ith, 4 ),
    ith_num( ith, 5 ) );
  fmt::print(
    "  {:<10} {:<20.15g} {:<20.15g} {:<20.15g} {:<20.15g} {:<20.15g} {:<20.15g}\n",
    "thetaM",
    thetaM( ith, 0 ),
    thetaM( ith, 1 ),
    thetaM( ith, 2 ),
    thetaM( ith, 3 ),
    thetaM( ith, 4 ),
    thetaM( ith, 5 ) );
  // average time
  fmt::print(
    "  {:<10} {:<20.15g} {:<20.15g} {:<20.15g} {:<20.15g} {:<20.15g} {:<20.15g}\n",
    "avg time",
    times.col( 0 ).mean(),
    times.col( 1 ).mean(),
    times.col( 2 ).mean(),
    times.col( 3 ).mean(),
    times.col( 4 ).mean(),
    times.col( 5 ).mean() );
  // average eval
  fmt::print(
    "  {:<10} {:<20.15g} {:<20.15g} {:<20.15g} {:<20.15g} {:<20.15g} {:<20.15g}\n",
    "avg eval",
    ith_num.col( 0 ).mean(),
    ith_num.col( 1 ).mean(),
    ith_num.col( 2 ).mean(),
    ith_num.col( 3 ).mean(),
    ith_num.col( 4 ).mean(),
    ith_num.col( 5 ).mean() );
  fmt::print( "---------------------------------------------\n" );
}

void printStatsFinal(
  Eigen::MatrixXd & times,
  Eigen::MatrixXd & ith_num,
  Eigen::MatrixXd & /* thetaM */
)
{
  fmt::print( "---------------------------------------------\n" );
  fmt::print(
    "  {:<10} {:<20} {:<20} {:<20} {:<20} {:<20} {:<20}\n",
    "VAL",
    "SD ",
    "PS ",
    "PS_BRACKET ",
    "PT ",
    "PT_BRACKET ",
    "ELL" );
  fmt::print( "  {:<10} {:<20} {:<20} {:<20} {:<20} {:<20} {:<20}\n", "---", "---", "---", "---", "---", "---", "---" );
  // average time
  fmt::print(
    "  {:<10} {:<20.15g} {:<20.15g} {:<20.15g} {:<20.15g} {:<20.15g} {:<20.15g}\n",
    "avg time",
    times.col( 0 ).mean(),
    times.col( 1 ).mean(),
    times.col( 2 ).mean(),
    times.col( 3 ).mean(),
    times.col( 4 ).mean(),
    times.col( 5 ).mean() );
  // average eval
  fmt::print(
    "  {:<10} {:<20.15g} {:<20.15g} {:<20.15g} {:<20.15g} {:<20.15g} {:<20.15g}\n",
    "avg eval",
    ith_num.col( 0 ).mean(),
    ith_num.col( 1 ).mean(),
    ith_num.col( 2 ).mean(),
    ith_num.col( 3 ).mean(),
    ith_num.col( 4 ).mean(),
    ith_num.col( 5 ).mean() );

  fmt::print(
    "  {:<10} {:<20.15g} {:<20.15g} {:<20.15g} {:<20.15g} {:<20.15g} {:<20.15g}\n",
    "t wrt 360",
    times.col( 0 ).mean() / times.col( 0 ).mean(),
    times.col( 0 ).mean() / times.col( 1 ).mean(),
    times.col( 0 ).mean() / times.col( 2 ).mean(),
    times.col( 0 ).mean() / times.col( 3 ).mean(),
    times.col( 0 ).mean() / times.col( 4 ).mean(),
    times.col( 0 ).mean() / times.col( 5 ).mean() );

  fmt::print(
    "  {:<10} {:<20.15g} {:<20.15g} {:<20.15g} {:<20.15g} {:<20.15g} {:<20.15g}\n",
    "it wrt 360",
    ith_num.col( 0 ).mean() / ith_num.col( 0 ).mean(),
    ith_num.col( 0 ).mean() / ith_num.col( 1 ).mean(),
    ith_num.col( 0 ).mean() / ith_num.col( 2 ).mean(),
    ith_num.col( 0 ).mean() / ith_num.col( 3 ).mean(),
    ith_num.col( 0 ).mean() / ith_num.col( 4 ).mean(),
    ith_num.col( 0 ).mean() / ith_num.col( 5 ).mean() );
  fmt::print( "---------------------------------------------\n" );
}
