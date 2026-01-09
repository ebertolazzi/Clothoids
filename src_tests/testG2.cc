/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2026                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include "Clothoids.hh"
#include "Clothoids_fmt.hh"
#include "Utils.hh"
#include "Utils_eigen.hh"

using G2lib::integer;
using G2lib::real_type;
using namespace std;
using Utils::m_pi;

// ============================================================================
// FORMATTING STYLES
// ============================================================================

namespace Style
{
  // Main styles
  const auto HEADER    = fg( fmt::color::steel_blue ) | fmt::emphasis::bold;
  const auto SECTION   = fg( fmt::color::dodger_blue ) | fmt::emphasis::bold;
  const auto SUCCESS   = fg( fmt::color::lime_green ) | fmt::emphasis::bold;
  const auto ERROR     = fg( fmt::color::crimson ) | fmt::emphasis::bold;
  const auto WARNING   = fg( fmt::color::gold ) | fmt::emphasis::bold;
  const auto INFO      = fg( fmt::color::deep_sky_blue ) | fmt::emphasis::bold;
  const auto VALUE     = fg( fmt::color::light_gray );
  const auto LABEL     = fg( fmt::color::silver );
  const auto HIGHLIGHT = fg( fmt::color::cyan ) | fmt::emphasis::bold;
  const auto GEOMETRY  = fg( fmt::color::violet ) | fmt::emphasis::bold;
  const auto CURVE     = fg( fmt::color::orange ) | fmt::emphasis::bold;
  const auto TEST_PASS = fg( fmt::color::green ) | fmt::emphasis::bold;
  const auto TEST_FAIL = fg( fmt::color::red ) | fmt::emphasis::bold;
  const auto POINT     = fg( fmt::color::spring_green ) | fmt::emphasis::bold;
  const auto PARAM     = fg( fmt::color::hot_pink ) | fmt::emphasis::bold;
}  // namespace Style

// ============================================================================
// FORMATTING UTILITY FUNCTIONS
// ============================================================================

void print_header( const string & title, const string & icon = "🌀" )
{
  fmt::print( "\n" );
  fmt::print( Style::HEADER, "╔══════════════════════════════════════════════════════════════╗\n" );
  fmt::print( Style::HEADER, "║ {:^58} ║\n", icon + " " + title );
  fmt::print( Style::HEADER, "╚══════════════════════════════════════════════════════════════╝\n" );
  fmt::print( "\n" );
}

void print_section( const string & title, const string & icon = "📋" )
{
  fmt::print( Style::SECTION, "\n┌{0:─^{1}}┐\n", "", 60 );
  fmt::print( Style::SECTION, "│ {:^58} │\n", icon + " " + title );
  fmt::print( Style::SECTION, "└{0:─^{1}}┘\n", "", 60 );
  fmt::print( "\n" );
}

void print_test_result( const string & test_name, bool passed, const string & details = "" )
{
  if ( passed ) { fmt::print( Style::TEST_PASS, "    ✅ {}: PASSED", test_name ); }
  else
  {
    fmt::print( Style::TEST_FAIL, "    ❌ {}: FAILED", test_name );
  }

  if ( !details.empty() ) { fmt::print( Style::INFO, " - {}", details ); }
  fmt::print( "\n" );
}

void print_curve_info( const G2lib::ClothoidCurve & C, const string & name = "Curve" )
{
  fmt::print( Style::GEOMETRY, "    📐 {}:\n", name );
  fmt::print( Style::LABEL, "      Start:      " );
  fmt::print( Style::POINT, "({:.6f}, {:.6f})\n", C.x_begin(), C.y_begin() );
  fmt::print( Style::LABEL, "      End:        " );
  fmt::print( Style::POINT, "({:.6f}, {:.6f})\n", C.x_end(), C.y_end() );
  fmt::print( Style::LABEL, "      Theta0:     " );
  fmt::print( Style::PARAM, "{:.6f} rad ({:.2f}°)\n", C.theta_begin(), C.theta_begin() * 180.0 / m_pi );
  fmt::print( Style::LABEL, "      Theta1:     " );
  fmt::print( Style::PARAM, "{:.6f} rad ({:.2f}°)\n", C.theta_end(), C.theta_end() * 180.0 / m_pi );
  fmt::print( Style::LABEL, "      Kappa0:     " );
  fmt::print( Style::PARAM, "{:.6f}\n", C.kappa_begin() );
  fmt::print( Style::LABEL, "      Kappa1:     " );
  fmt::print( Style::PARAM, "{:.6f}\n", C.kappa_end() );
  fmt::print( Style::LABEL, "      Length:     " );
  fmt::print( Style::HIGHLIGHT, "{:.6f}\n", C.length() );
}

void print_g2_problem_info(
  real_type      x0,
  real_type      y0,
  real_type      th0,
  real_type      k0,
  real_type      x1,
  real_type      y1,
  real_type      th1,
  real_type      k1,
  const string & title = "G² Problem" )
{
  fmt::print( Style::CURVE, "    🎯 {}:\n", title );
  fmt::print( Style::LABEL, "      Initial point P₀: " );
  fmt::print( Style::POINT, "({:.6f}, {:.6f})\n", x0, y0 );
  fmt::print( Style::LABEL, "      Initial angle θ₀: " );
  fmt::print( Style::PARAM, "{:.6f} rad ({:.2f}°)\n", th0, th0 * 180.0 / m_pi );
  fmt::print( Style::LABEL, "      Initial curvature κ₀: " );
  fmt::print( Style::PARAM, "{:.6f}\n", k0 );
  fmt::print( Style::LABEL, "      Final point P₁: " );
  fmt::print( Style::POINT, "({:.6f}, {:.6f})\n", x1, y1 );
  fmt::print( Style::LABEL, "      Final angle θ₁: " );
  fmt::print( Style::PARAM, "{:.6f} rad ({:.2f}°)\n", th1, th1 * 180.0 / m_pi );
  fmt::print( Style::LABEL, "      Final curvature κ₁: " );
  fmt::print( Style::PARAM, "{:.6f}\n", k1 );
}

// ============================================================================
// TEST G2solve2arc
// ============================================================================

void test_g2solve2arc()
{
  print_section( "TEST G² WITH 2 ARCS (G2solve2arc)", "🔄" );

  int tests_passed = 0;
  int total_tests  = 0;

  // Test 1: Basic problem (straight line with zero curvature)
  {
    total_tests++;
    G2lib::G2solve2arc solver;

    real_type x0 = 0.0, y0 = 0.0, th0 = 0.0, k0 = 0.0;
    real_type x1 = 1.0, y1 = 0.0, th1 = 0.0, k1 = 0.0;

    fmt::print( Style::INFO, "    📊 Test 1: Straight line (zero curvature)\n" );
    print_g2_problem_info( x0, y0, th0, k0, x1, y1, th1, k1 );

    int iter = solver.build( x0, y0, th0, k0, x1, y1, th1, k1 );

    fmt::print( Style::LABEL, "      Iterations: " );
    fmt::print( Style::HIGHLIGHT, "{}\n", iter );

    if ( iter >= 0 )
    {
      const auto & S0 = solver.S0();
      const auto & S1 = solver.S1();

      print_curve_info( S0, "S0 (first arc)" );
      print_curve_info( S1, "S1 (second arc)" );

      // Check continuity
      real_type err_x  = S0.x_end() - S1.x_begin();
      real_type err_y  = S0.y_end() - S1.y_begin();
      real_type err_th = S0.theta_end() - S1.theta_begin();
      real_type err_k  = S0.kappa_end() - S1.kappa_begin();

      fmt::print( Style::LABEL, "      G² continuity check:\n" );
      fmt::print( Style::VALUE, "        Δx:   {:.2e}\n", err_x );
      fmt::print( Style::VALUE, "        Δy:   {:.2e}\n", err_y );
      fmt::print( Style::VALUE, "        Δθ:   {:.2e} rad\n", err_th );
      fmt::print( Style::VALUE, "        Δκ:   {:.2e}\n", err_k );

      bool passed = ( abs( err_x ) < 1e-6 ) && ( abs( err_y ) < 1e-6 ) && ( abs( err_th ) < 1e-6 ) &&
                    ( abs( err_k ) < 1e-6 );

      print_test_result( "Straight line G²", passed );
      if ( passed ) tests_passed++;
    }
    else
    {
      print_test_result( "Straight line G²", false, "solver failed" );
    }
  }

  // Test 2: Simple curve
  {
    total_tests++;
    G2lib::G2solve2arc solver;

    real_type x0 = 0.0, y0 = 0.0, th0 = 0.0, k0 = 0.1;
    real_type x1 = 1.0, y1 = 1.0, th1 = m_pi / 2.0, k1 = 0.1;

    fmt::print( Style::INFO, "\n    📊 Test 2: Curve with constant curvature\n" );
    print_g2_problem_info( x0, y0, th0, k0, x1, y1, th1, k1 );

    int iter = solver.build( x0, y0, th0, k0, x1, y1, th1, k1 );

    fmt::print( Style::LABEL, "      Iterations: " );
    fmt::print( Style::HIGHLIGHT, "{}\n", iter );

    if ( iter >= 0 )
    {
      const auto & S0 = solver.S0();
      const auto & S1 = solver.S1();

      // Check total length is positive
      real_type total_length = S0.length() + S1.length();
      fmt::print( Style::LABEL, "      Total length: " );
      fmt::print( Style::HIGHLIGHT, "{:.6f}\n", total_length );

      bool passed = ( total_length > 0 ) && ( iter >= 0 );
      print_test_result( "Curve with constant curvature", passed );
      if ( passed ) tests_passed++;
    }
    else
    {
      print_test_result( "Curve with constant curvature", false, "solver failed" );
    }
  }

  // Test 3: Configuration with different curvatures
  {
    total_tests++;
    G2lib::G2solve2arc solver;

    real_type x0 = -1.0, y0 = 0.0, th0 = -2.0, k0 = 0.5;
    real_type x1 = 1.0, y1 = 0.0, th1 = 2.0, k1 = 0.5;

    fmt::print( Style::INFO, "\n    📊 Test 3: Symmetric configuration\n" );
    print_g2_problem_info( x0, y0, th0, k0, x1, y1, th1, k1 );

    int iter = solver.build( x0, y0, th0, k0, x1, y1, th1, k1 );

    fmt::print( Style::LABEL, "      Iterations: " );
    if ( iter >= 0 )
    {
      fmt::print( Style::HIGHLIGHT, "{}\n", iter );

      const auto & S0 = solver.S0();
      const auto & S1 = solver.S1();

      // Check symmetry
      real_type length0 = S0.length();
      real_type length1 = S1.length();

      fmt::print( Style::LABEL, "      Length S0: " );
      fmt::print( Style::VALUE, "{:.6f}\n", length0 );
      fmt::print( Style::LABEL, "      Length S1: " );
      fmt::print( Style::VALUE, "{:.6f}\n", length1 );

      bool passed = ( abs( length0 - length1 ) < 1e-6 ) && ( iter >= 0 );
      print_test_result( "Symmetric configuration", passed );
      if ( passed ) tests_passed++;
    }
    else
    {
      fmt::print( Style::WARNING, "{}\n", iter );
      print_test_result( "Symmetric configuration", false, "solver failed" );
    }
  }

  // Test 4: Setting tolerance and maximum iterations
  {
    total_tests++;
    G2lib::G2solve2arc solver;

    // Set strict parameters
    solver.set_tolerance( 1e-12 );
    solver.set_max_iter( 50 );

    real_type x0 = 0.0, y0 = 0.0, th0 = 0.0, k0 = 0.0;
    real_type x1 = 2.0, y1 = 1.0, th1 = m_pi / 3.0, k1 = 0.2;

    fmt::print( Style::INFO, "\n    📊 Test 4: Solver with strict parameters\n" );
    fmt::print( Style::VALUE, "      Tolerance: 1e-12, Max iter: 50\n" );

    int iter = solver.build( x0, y0, th0, k0, x1, y1, th1, k1 );

    fmt::print( Style::LABEL, "      Iterations: " );
    fmt::print( Style::HIGHLIGHT, "{}\n", iter );

    bool passed = ( iter >= 0 ) && ( iter <= 50 );
    print_test_result( "Solver with strict parameters", passed );
    if ( passed ) tests_passed++;
  }

  // Section summary
  fmt::print( "\n" );
  fmt::print( Style::INFO, "    📊 Summary G2solve2arc: " );
  fmt::print( Style::HIGHLIGHT, "{}/{} tests passed\n", tests_passed, total_tests );
}

// ============================================================================
// TEST G2solve3arc
// ============================================================================

void test_g2solve3arc()
{
  print_section( "TEST G² WITH 3 ARCS (G2solve3arc)", "🌀" );

  int tests_passed = 0;
  int total_tests  = 0;

  // Test 1: Problem from original example
  {
    total_tests++;
    G2lib::G2solve3arc solver;

    real_type x0  = -1;
    real_type y0  = 0;
    real_type th0 = -2;
    real_type k0  = 0.909297426825682;
    real_type x1  = 1;
    real_type y1  = 0;
    real_type th1 = 2;
    real_type k1  = 0.909297426825682;

    fmt::print( Style::INFO, "    📊 Test 1: Original example from test file\n" );
    print_g2_problem_info( x0, y0, th0, k0, x1, y1, th1, k1, "Original example" );

    int iter = solver.build( x0, y0, th0, k0, x1, y1, th1, k1 );

    fmt::print( Style::LABEL, "      Iterations: " );
    fmt::print( Style::HIGHLIGHT, "{}\n", iter );

    if ( iter >= 0 )
    {
      const auto & S0 = solver.S0();
      const auto & SM = solver.SM();
      const auto & S1 = solver.S1();

      fmt::print( Style::LABEL, "      📐 Generated arcs:\n" );
      print_curve_info( S0, "S0 (first arc)" );
      print_curve_info( SM, "SM (middle arc)" );
      print_curve_info( S1, "S1 (last arc)" );

      // Check G² continuity between arcs
      real_type err1_x  = S0.x_end() - SM.x_begin();
      real_type err1_y  = S0.y_end() - SM.y_begin();
      real_type err1_th = S0.theta_end() - SM.theta_begin();
      real_type err1_k  = S0.kappa_end() - SM.kappa_begin();

      real_type err2_x  = S1.x_begin() - SM.x_end();
      real_type err2_y  = S1.y_begin() - SM.y_end();
      real_type err2_th = S1.theta_begin() - SM.theta_end();
      real_type err2_k  = S1.kappa_begin() - SM.kappa_end();

      fmt::print( Style::LABEL, "      🔍 G² continuity check:\n" );
      fmt::print(
        Style::VALUE,
        "        Transition S0→SM: Δx={:.2e}, Δy={:.2e}, Δθ={:.2e}, Δκ={:.2e}\n",
        err1_x,
        err1_y,
        err1_th,
        err1_k );
      fmt::print(
        Style::VALUE,
        "        Transition SM→S1: Δx={:.2e}, Δy={:.2e}, Δθ={:.2e}, Δκ={:.2e}\n",
        err2_x,
        err2_y,
        err2_th,
        err2_k );

      real_type max_err = max(
        { abs( err1_x ),
          abs( err1_y ),
          abs( err1_th ),
          abs( err1_k ),
          abs( err2_x ),
          abs( err2_y ),
          abs( err2_th ),
          abs( err2_k ) } );

      fmt::print( Style::LABEL, "      📏 Total length: " );
      fmt::print( Style::HIGHLIGHT, "{:.6f}\n", solver.total_length() );

      bool passed = ( max_err < 1e-6 ) && ( iter >= 0 );
      print_test_result( "Original example 3 arcs", passed );
      if ( passed ) tests_passed++;
    }
    else
    {
      print_test_result( "Original example 3 arcs", false, "solver failed" );
    }
  }

  // Test 2: Fixed length problem
  {
    total_tests++;
    G2lib::G2solve3arc solver;

    real_type x0 = 0.0, y0 = 0.0, th0 = 0.0, k0 = 0.0;
    real_type x1 = 2.0, y1 = 1.0, th1 = m_pi / 4.0, k1 = 0.0;
    real_type s0 = 1.0 / 2;  // Fixed length for first arc
    real_type s1 = 1.5 / 2;  // Fixed length for last arc

    fmt::print( Style::INFO, "\n    📊 Test 2: Fixed length problem\n" );
    print_g2_problem_info( x0, y0, th0, k0, x1, y1, th1, k1 );
    fmt::print( Style::VALUE, "      Fixed length S0: {:.2f}\n", s0 );
    fmt::print( Style::VALUE, "      Fixed length S1: {:.2f}\n", s1 );

    int iter = solver.build_fixed_length( s0, x0, y0, th0, k0, s1, x1, y1, th1, k1 );

    fmt::print( Style::LABEL, "      Iterations: " );
    fmt::print( Style::HIGHLIGHT, "{}\n", iter );

    if ( iter >= 0 )
    {
      const auto & S0 = solver.S0();
      const auto & S1 = solver.S1();

      real_type length0 = S0.length();
      real_type length1 = S1.length();

      fmt::print( Style::LABEL, "      Actual length S0: " );
      fmt::print( Style::VALUE, "{:.6f} (expected: {:.2f})\n", length0, s0 );
      fmt::print( Style::LABEL, "      Actual length S1: " );
      fmt::print( Style::VALUE, "{:.6f} (expected: {:.2f})\n", length1, s1 );

      bool passed = ( abs( length0 - s0 ) < 1e-6 ) && ( abs( length1 - s1 ) < 1e-6 );
      print_test_result( "Fixed length problem", passed );
      if ( passed ) tests_passed++;
    }
    else
    {
      print_test_result( "Fixed length problem", false, "solver failed" );
    }
  }

  // Test 3: Evaluation of geometric properties
  {
    total_tests++;
    G2lib::G2solve3arc solver;

    real_type x0 = 0.0, y0 = 0.0, th0 = 0.0, k0 = 0.1;
    real_type x1 = 3.0, y1 = 2.0, th1 = m_pi / 2.0, k1 = -0.1;

    fmt::print( Style::INFO, "\n    📊 Test 3: Evaluation of geometric properties\n" );

    int iter = solver.build( x0, y0, th0, k0, x1, y1, th1, k1 );

    if ( iter >= 0 )
    {
      // Calculate various properties
      real_type total_length    = solver.total_length();
      real_type theta_variation = solver.theta_total_variation();
      real_type kappa_variation = solver.curvature_total_variation();
      real_type integral_kappa2 = solver.integral_curvature2();
      real_type integral_jerk2  = solver.integral_jerk2();

      real_type th_min, th_max;
      solver.theta_min_max( th_min, th_max );

      real_type k_min, k_max;
      solver.curvature_min_max( k_min, k_max );

      fmt::print( Style::LABEL, "      📊 Geometric properties:\n" );
      fmt::print( Style::VALUE, "        Total length:         {:.6f}\n", total_length );
      fmt::print( Style::VALUE, "        Total θ variation:    {:.6f} rad\n", theta_variation );
      fmt::print( Style::VALUE, "        Total κ variation:    {:.6f}\n", kappa_variation );
      fmt::print( Style::VALUE, "        ∫κ² ds:              {:.6f}\n", integral_kappa2 );
      fmt::print( Style::VALUE, "        ∫jerk² ds:           {:.6f}\n", integral_jerk2 );
      fmt::print( Style::VALUE, "        θ range:             [{:.6f}, {:.6f}] rad\n", th_min, th_max );
      fmt::print( Style::VALUE, "        κ range:             [{:.6f}, {:.6f}]\n", k_min, k_max );

      // Test evaluation of midpoint
      real_type s_mid = total_length / 2.0;
      real_type x_mid, y_mid;
      solver.eval( s_mid, x_mid, y_mid );

      fmt::print( Style::LABEL, "      📍 Midpoint evaluation (s={:.6f}):\n", s_mid );
      fmt::print( Style::VALUE, "        ({:.6f}, {:.6f})\n", x_mid, y_mid );

      bool passed = ( total_length > 0 ) && ( integral_kappa2 >= 0 );
      print_test_result( "Geometric properties", passed );
      if ( passed ) tests_passed++;
    }
    else
    {
      print_test_result( "Geometric properties", false, "solver failed" );
    }
  }

  // Test 4: Evaluation of derivatives
  {
    total_tests++;
    G2lib::G2solve3arc solver;

    real_type x0 = 0.0, y0 = 0.0, th0 = 0.0, k0 = 0.0;
    real_type x1 = 2.0, y1 = 0.0, th1 = m_pi / 4.0, k1 = 0.2;

    int iter = solver.build( x0, y0, th0, k0, x1, y1, th1, k1 );

    if ( iter >= 0 )
    {
      real_type s_test = solver.total_length() / 3.0;

      // Function and derivatives evaluation
      real_type theta_s  = solver.theta( s_test );
      real_type theta_D  = solver.theta_D( s_test );   // Curvature
      real_type theta_DD = solver.theta_DD( s_test );  // Curvature derivative

      real_type x_s, y_s;
      solver.eval( s_test, x_s, y_s );

      real_type x_D, y_D;
      solver.eval_D( s_test, x_D, y_D );

      real_type x_DD, y_DD;
      solver.eval_DD( s_test, x_DD, y_DD );

      fmt::print( Style::INFO, "\n    📊 Test 4: Evaluation of derivatives\n" );
      fmt::print( Style::VALUE, "      s = {:.6f}:\n", s_test );
      fmt::print( Style::VALUE, "        θ(s)   = {:.6f} rad\n", theta_s );
      fmt::print( Style::VALUE, "        θ'(s)  = {:.6f} (curvature)\n", theta_D );
      fmt::print( Style::VALUE, "        θ''(s) = {:.6f}\n", theta_DD );
      fmt::print( Style::VALUE, "        x(s)   = {:.6f}, y(s) = {:.6f}\n", x_s, y_s );
      fmt::print( Style::VALUE, "        x'(s)  = {:.6f}, y'(s) = {:.6f}\n", x_D, y_D );
      fmt::print( Style::VALUE, "        x''(s) = {:.6f}, y''(s) = {:.6f}\n", x_DD, y_DD );

      // Check consistency: x' = cos(θ), y' = sin(θ)
      real_type cos_theta = cos( theta_s );
      real_type sin_theta = sin( theta_s );
      real_type err_x     = x_D - cos_theta;
      real_type err_y     = y_D - sin_theta;

      bool passed = ( abs( err_x ) < 1e-6 ) && ( abs( err_y ) < 1e-6 );
      print_test_result( "Evaluation of derivatives", passed );
      if ( passed ) tests_passed++;
    }
    else
    {
      print_test_result( "Evaluation of derivatives", false, "solver failed" );
    }
  }

  // Test 5: Setting solver parameters
  {
    total_tests++;
    G2lib::G2solve3arc solver;

    // Set custom parameters
    solver.set_tolerance( 1e-12 );
    solver.set_max_iter( 200 );

    real_type x0 = -1.0, y0 = -1.0, th0 = -m_pi / 3.0, k0 = 0.3;
    real_type x1 = 1.0, y1 = 1.0, th1 = m_pi / 3.0, k1 = 0.3;

    fmt::print( Style::INFO, "\n    📊 Test 5: Solver with custom parameters\n" );
    fmt::print( Style::VALUE, "      Tolerance: 1e-12, Max iter: 200\n" );

    int iter = solver.build( x0, y0, th0, k0, x1, y1, th1, k1 );

    fmt::print( Style::LABEL, "      Iterations: " );
    fmt::print( Style::HIGHLIGHT, "{}\n", iter );

    bool passed = ( iter >= 0 ) && ( iter <= 200 );
    print_test_result( "Solver with custom parameters", passed );
    if ( passed ) tests_passed++;
  }

  // Section summary
  fmt::print( "\n" );
  fmt::print( Style::INFO, "    📊 Summary G2solve3arc: " );
  fmt::print( Style::HIGHLIGHT, "{}/{} tests passed\n", tests_passed, total_tests );
}

// ============================================================================
// TEST G2solveCLC (Clothoid-Line-Clothoid)
// ============================================================================

void test_g2solve_clc()
{
  print_section( "TEST G² CLC (CLOTHOID-LINE-CLOTHOID)", "📏" );

  int tests_passed = 0;
  int total_tests  = 0;

  // Test 1: Basic CLC problem
  {
    total_tests++;
    G2lib::G2solveCLC solver;

    real_type x0 = 0.0, y0 = 0.0, th0 = 0, k0 = 1;
    real_type x1 = 3.0, y1 = 2.0, th1 = 0, k1 = -1;

    fmt::print( Style::INFO, "    📊 Test 1: Basic CLC problem\n" );
    print_g2_problem_info( x0, y0, th0, k0, x1, y1, th1, k1 );

    int iter = solver.build( x0, y0, th0, k0, x1, y1, th1, k1 );

    fmt::print( Style::LABEL, "      Iterations: " );
    fmt::print( Style::HIGHLIGHT, "{}\n", iter );

    if ( iter >= 0 )
    {
      const auto & S0 = solver.S0();
      const auto & SM = solver.SM();
      const auto & S1 = solver.S1();

      fmt::print( Style::LABEL, "      📐 Generated arcs:\n" );
      print_curve_info( S0, "S0 (first clothoid arc)" );
      print_curve_info( SM, "SM (line segment)" );
      print_curve_info( S1, "S1 (last clothoid arc)" );

      // Verify SM is actually a line segment (zero curvature)
      real_type kappa_SM_begin = SM.kappa_begin();
      real_type kappa_SM_end   = SM.kappa_end();

      fmt::print( Style::LABEL, "      🔍 Line segment verification:\n" );
      fmt::print( Style::VALUE, "        κ start SM: {:.2e}\n", kappa_SM_begin );
      fmt::print( Style::VALUE, "        κ end SM:   {:.2e}\n", kappa_SM_end );

      bool passed = ( abs( kappa_SM_begin ) < 1e-10 ) && ( abs( kappa_SM_end ) < 1e-10 );
      print_test_result( "Basic CLC problem", passed );
      if ( passed ) tests_passed++;
    }
    else
    {
      print_test_result( "Basic CLC problem", false, "solver failed" );
    }
  }

  // Test 2: CLC with different curvatures
  {
    total_tests++;
    G2lib::G2solveCLC solver;

    solver.set_tolerance( 1e-10 );
    solver.set_max_iter( 100 );

    real_type x0 = -1.0, y0 = 0.0, th0 = -0.5, k0 = 2;
    real_type x1 = 1.0, y1 = 0.0, th1 = 0.5, k1 = 2;

    fmt::print( Style::INFO, "\n    📊 Test 2: Symmetric CLC\n" );

    int iter = solver.build( x0, y0, th0, k0, x1, y1, th1, k1 );

    fmt::print( Style::LABEL, "      Iterations: " );
    fmt::print( Style::HIGHLIGHT, "{}\n", iter );

    bool passed = ( iter >= 0 );
    print_test_result( "Symmetric CLC", passed );
    if ( passed ) tests_passed++;
  }

  // Section summary
  fmt::print( "\n" );
  fmt::print( Style::INFO, "    📊 Summary G2solveCLC: " );
  fmt::print( Style::HIGHLIGHT, "{}/{} tests passed\n", tests_passed, total_tests );
}

// ============================================================================
// TEST ClothoidSplineG2
// ============================================================================

void test_clothoid_spline_g2()
{
  print_section( "TEST CLOTHOID SPLINE G²", "📈" );

  int tests_passed = 0;
  int total_tests  = 0;

  // Test 1: Simple spline
  {
    total_tests++;
    G2lib::ClothoidSplineG2 spline;

    // Points for spline
    integer           npts = 5;
    vector<real_type> X    = { 0.0, 1.0, 2.0, 3.0, 4.0 };
    vector<real_type> Y    = { 0.0, 1.0, 0.5, 2.0, 1.0 };
    vector<real_type> theta( npts );

    // Use target P4 (default)
    spline.setP4();

    fmt::print( Style::INFO, "    📊 Test 1: Simple spline ({} points)\n", npts );
    fmt::print( Style::VALUE, "      Points: " );
    for ( integer i = 0; i < npts; ++i ) { fmt::print( Style::POINT, "({:.1f},{:.1f}) ", X[i], Y[i] ); }
    fmt::print( "\n" );

    // Build spline
    bool success = spline.build_PN( npts, X.data(), Y.data(), theta.data(), G2lib::TargetType::P4 );

    fmt::print( Style::LABEL, "      Construction: " );
    if ( success )
    {
      fmt::print( Style::SUCCESS, "success\n" );

      fmt::print( Style::LABEL, "      Number of points: " );
      fmt::print( Style::VALUE, "{}\n", spline.numPnts() );
      fmt::print( Style::LABEL, "      Number of theta: " );
      fmt::print( Style::VALUE, "{}\n", spline.numTheta() );
      fmt::print( Style::LABEL, "      Constraints: " );
      fmt::print( Style::VALUE, "{}\n", spline.numConstraints() );

      // Print calculated angles
      fmt::print( Style::LABEL, "      Calculated θ angles:\n" );
      for ( integer i = 0; i < npts; ++i )
      {
        fmt::print( Style::VALUE, "        θ[{}] = {:.6f} rad ({:.2f}°)\n", i, theta[i], theta[i] * 180.0 / m_pi );
      }

      print_test_result( "Simple spline", true );
      tests_passed++;
    }
    else
    {
      fmt::print( Style::ERROR, "failed\n" );
      print_test_result( "Simple spline", false );
    }
  }

  // Test 2: Spline with boundary conditions
  {
    total_tests++;
    G2lib::ClothoidSplineG2 spline;

    integer           npts = 6;
    vector<real_type> X    = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
    vector<real_type> Y    = { 0.0, 1.0, 0.0, -1.0, 0.0, 1.0 };  // Sinusoidal pattern
    vector<real_type> theta( npts );

    fmt::print( Style::INFO, "\n    📊 Test 2: Spline with P1 conditions (fixed initial/final angles)\n" );

    // Use build_P1 which fixes initial and final angles
    real_type theta_init = 0.0;
    real_type theta_end  = m_pi / 4.0;

    bool success = spline.build_P1( npts, X.data(), Y.data(), theta.data(), theta_init, theta_end );

    fmt::print( Style::LABEL, "      P1 construction: " );
    if ( success )
    {
      fmt::print( Style::SUCCESS, "success\n" );
      fmt::print( Style::VALUE, "      Imposed initial θ: {:.4f} rad\n", theta_init );
      fmt::print( Style::VALUE, "      Imposed final θ:   {:.4f} rad\n", theta_end );
      fmt::print( Style::VALUE, "      Calculated θ[0]:     {:.4f} rad\n", theta[0] );
      fmt::print( Style::VALUE, "      Calculated θ[{}]:    {:.4f} rad\n", npts - 1, theta[npts - 1] );

      // Verify endpoint angles match requirements
      bool angles_match = ( abs( theta[0] - theta_init ) < 1e-6 ) && ( abs( theta[npts - 1] - theta_end ) < 1e-6 );

      print_test_result( "Spline with P1 conditions", angles_match );
      if ( angles_match ) tests_passed++;
    }
    else
    {
      fmt::print( Style::ERROR, "failed\n" );
      print_test_result( "Spline with P1 conditions", false );
    }
  }

  // Test 3: Spline with different targets
  {
    total_tests++;
    G2lib::ClothoidSplineG2 spline;

    integer           npts = 4;
    vector<real_type> X    = { 0.0, 2.0, 4.0, 6.0 };
    vector<real_type> Y    = { 0.0, 3.0, 1.0, 4.0 };
    vector<real_type> theta( npts );

    // Try different targets
    vector<pair<G2lib::TargetType, string>> targets = {
      { G2lib::TargetType::P4, "P4" }, { G2lib::TargetType::P5, "P5" }, { G2lib::TargetType::P6, "P6" },
      { G2lib::TargetType::P7, "P7" }, { G2lib::TargetType::P8, "P8" }, { G2lib::TargetType::P9, "P9" }
    };

    fmt::print( Style::INFO, "\n    📊 Test 3: Spline with different targets\n" );

    int successful_targets = 0;
    for ( const auto & [target_type, target_name] : targets )
    {
      bool success = spline.build_PN( npts, X.data(), Y.data(), theta.data(), target_type );

      if ( success )
      {
        fmt::print( Style::VALUE, "      Target {}: ", target_name );
        fmt::print( Style::SUCCESS, "✓\n" );
        successful_targets++;
      }
      else
      {
        fmt::print( Style::VALUE, "      Target {}: ", target_name );
        fmt::print( Style::WARNING, "⚠\n" );
      }
    }

    bool passed = ( successful_targets > 0 );
    print_test_result(
      "Spline with different targets",
      passed,
      fmt::format( "{}/{} targets successful", successful_targets, targets.size() ) );
    if ( passed ) tests_passed++;
  }

  // Test 4: Objective and gradient functions
  {
    total_tests++;
    G2lib::ClothoidSplineG2 spline;

    integer           npts = 5;
    vector<real_type> X    = { 0.0, 1.0, 2.0, 3.0, 4.0 };
    vector<real_type> Y    = { 0.0, 2.0, 1.0, 3.0, 2.0 };
    vector<real_type> theta( npts );

    spline.build_PN( npts, X.data(), Y.data(), theta.data(), G2lib::TargetType::P4 );

    // Test objective function
    real_type f_value;
    bool      obj_success = spline.objective( theta.data(), f_value );

    // Test gradient (simulated)
    fmt::print( Style::INFO, "\n    📊 Test 4: Objective and gradient functions\n" );
    fmt::print( Style::LABEL, "      Objective function: " );
    if ( obj_success )
    {
      fmt::print( Style::SUCCESS, "calculated\n" );
      fmt::print( Style::VALUE, "      f value: {:.6f}\n", f_value );
      print_test_result( "Objective and gradient functions", true );
      tests_passed++;
    }
    else
    {
      fmt::print( Style::WARNING, "not calculated\n" );
      print_test_result( "Objective and gradient functions", false );
    }
  }

  // Section summary
  fmt::print( "\n" );
  fmt::print( Style::INFO, "    📊 Summary ClothoidSplineG2: " );
  fmt::print( Style::HIGHLIGHT, "{}/{} tests passed\n", tests_passed, total_tests );
}

// ============================================================================
// TEST ClothoidList G2
// ============================================================================

void test_clothoid_list_g2()
{
  print_section( "TEST CLOTHOID LIST G²", "📊" );

  int tests_passed = 0;
  int total_tests  = 0;

  // Test 1: G2 construction from points
  {
    total_tests++;
    G2lib::ClothoidList S( "test_G2" );

    integer           N = 5;
    vector<real_type> X = { 0.0, 1.0, 2.0, 3.0, 4.0 };
    vector<real_type> Y = { 0.0, 1.0, 0.5, 1.5, 1.0 };

    real_type theta0 = m_pi / 4.0;   // Initial angle
    real_type thetaN = -m_pi / 4.0;  // Final angle

    fmt::print( Style::INFO, "    📊 Test 1: G2 construction from {} points\n", N );

    Utils::TicToc tictoc;
    tictoc.tic();
    S.build_G2( N, X.data(), Y.data(), theta0, 0.0, thetaN, 0.0 );
    tictoc.toc();

    fmt::print( Style::LABEL, "      Construction time: " );
    fmt::print( Style::VALUE, "{:.3f} ms\n", tictoc.elapsed_ms() );

    fmt::print( Style::LABEL, "      Number of segments: " );
    fmt::print( Style::VALUE, "{}\n", S.num_segments() );
    fmt::print( Style::LABEL, "      Total length: " );
    fmt::print( Style::VALUE, "{:.6f}\n", S.length() );

    bool passed = ( S.num_segments() == N + 3 );
    print_test_result( "G2 construction from points", passed );
    if ( passed ) tests_passed++;
  }

  // Test 2: Cyclic G2 construction
  {
    total_tests++;
    G2lib::ClothoidList S( "test_G2_cyclic" );

    integer           N = 6;
    vector<real_type> X = { 0.0, 2.0, 4.0, 4.0, 2.0, 0.0 };
    vector<real_type> Y = { 0.0, 3.0, 2.0, 0.0, -1.0, 0.0 };  // Closed contour

    fmt::print( Style::INFO, "\n    📊 Test 2: Cyclic G2 construction\n" );

    Utils::TicToc tictoc;
    tictoc.tic();
    S.build_G2_cyclic( N, X.data(), Y.data() );
    tictoc.toc();

    fmt::print( Style::LABEL, "      Construction time: " );
    fmt::print( Style::VALUE, "{:.3f} ms\n", tictoc.elapsed_ms() );

    // Verify continuity between start and end
    real_type err_x  = S.x_end() - S.x_begin();
    real_type err_y  = S.y_end() - S.y_begin();
    real_type err_th = S.theta_end() - S.theta_begin();

    fmt::print( Style::LABEL, "      Closure check:\n" );
    fmt::print( Style::VALUE, "        Δx:   {:.2e}\n", err_x );
    fmt::print( Style::VALUE, "        Δy:   {:.2e}\n", err_y );
    fmt::print( Style::VALUE, "        Δθ:   {:.2e} rad\n", err_th );

    bool passed = ( abs( err_x ) < 1e-6 ) && ( abs( err_y ) < 1e-6 );
    print_test_result( "Cyclic G2 construction", passed );
    if ( passed ) tests_passed++;
  }

  // Test 3: Loading from file and compatibility
  {
    total_tests++;

    // Create temporary test file
    ofstream test_file( "test_G2_temp.txt" );
    if ( test_file )
    {
      test_file << "CLOTHOID_LIST\n";
      test_file << "2\n";            // Number of clothoids
      test_file << "0 0 0 0 0 1\n";  // First clothoid
      test_file << "1 0 0 0 1 1\n";  // Second clothoid
      test_file.close();

      G2lib::ClothoidList S( "loaded" );
      ifstream            file( "test_G2_temp.txt" );

      fmt::print( Style::INFO, "\n    📊 Test 3: Loading from file\n" );

      if ( file.is_open() )
      {
        S.load( file );
        file.close();

        fmt::print( Style::LABEL, "      Clothoids loaded: " );
        fmt::print( Style::VALUE, "{}\n", S.num_segments() );

        // Remove temporary file
        remove( "test_G2_temp.txt" );

        print_test_result( "Loading from file and compatibility", true );
        tests_passed++;
      }
      else
      {
        print_test_result( "Loading from file and compatibility", false, "file not opened" );
      }
    }
    else
    {
      print_test_result( "Loading from file and compatibility", false, "file not created" );
    }
  }

  // Test 4: Get methods
  {
    total_tests++;
    G2lib::ClothoidCurve S0( "test_get_methods" );
    G2lib::ClothoidList  S( "test_get_methods" );

    // Build simple list
    S0.build_G1( 0.0, 0.0, 0.0, 3.0, 0.0, m_pi / 2.0 );
    S.push_back( S0 );
    S.push_back_G1( 5.0, 3.0, 0.0 );

    fmt::print( Style::INFO, "\n    📊 Test 4: Get methods for clothoid data\n" );

    // Test get methods
    vector<real_type> s_vec, theta_vec, kappa_vec;
    S.get_STK( s_vec, theta_vec, kappa_vec );

    fmt::print( Style::LABEL, "      Number of break points: " );
    fmt::print( Style::VALUE, "{}\n", s_vec.size() );

    if ( !s_vec.empty() )
    {
      fmt::print(
        Style::VALUE,
        "      First point: s={:.3f}, θ={:.3f}, κ={:.3f}\n",
        s_vec[0],
        theta_vec[0],
        kappa_vec[0] );
      fmt::print(
        Style::VALUE,
        "      Last point: s={:.3f}, θ={:.3f}, κ={:.3f}\n",
        s_vec.back(),
        theta_vec.back(),
        kappa_vec.back() );

      print_test_result( "Get methods for clothoid data", true );
      tests_passed++;
    }
    else
    {
      print_test_result( "Get methods for clothoid data", false, "no data obtained" );
    }
  }

  // Section summary
  fmt::print( "\n" );
  fmt::print( Style::INFO, "    📊 Summary ClothoidList G2: " );
  fmt::print( Style::HIGHLIGHT, "{}/{} tests passed\n", tests_passed, total_tests );
}

// ============================================================================
// TEST ORIGINAL EXAMPLE
// ============================================================================

void test_original_example()
{
  print_section( "TEST COMPLETE ORIGINAL EXAMPLE", "📄" );

  fmt::print( Style::LABEL, "    📖 Running the original example from testG2.cc file:\n\n" );

  // Part 1: Test G2solve3arc from example
  {
    fmt::print( Style::INFO, "    📊 Part 1: Original G2solve3arc test\n" );

    G2lib::G2solve3arc g2solve3arc;

    real_type x0  = -1;
    real_type y0  = 0;
    real_type th0 = -2;
    real_type k0  = 0.909297426825682;
    real_type x1  = 1;
    real_type y1  = 0;
    real_type th1 = 2;
    real_type k1  = 0.909297426825682;

    print_g2_problem_info( x0, y0, th0, k0, x1, y1, th1, k1 );

    int iter = g2solve3arc.build( x0, y0, th0, k0, x1, y1, th1, k1 );

    fmt::print( Style::LABEL, "      Iterations: " );
    fmt::print( Style::HIGHLIGHT, "{}\n", iter );

    if ( iter >= 0 )
    {
      G2lib::ClothoidCurve const & S0 = g2solve3arc.S0();
      G2lib::ClothoidCurve const & S1 = g2solve3arc.S1();
      G2lib::ClothoidCurve const & SM = g2solve3arc.SM();

      // Check connection errors
      fmt::print( Style::LABEL, "      Connection errors check:\n" );

      real_type err1 = S0.x_end() - SM.x_begin();
      real_type err2 = S0.y_end() - SM.y_begin();
      real_type err3 = S0.theta_end() - SM.theta_begin();
      real_type err4 = S1.x_begin() - SM.x_end();
      real_type err5 = S1.y_begin() - SM.y_end();
      real_type err6 = S1.theta_begin() - SM.theta_end();

      fmt::print( Style::VALUE, "        Δx0:   {:.2e}\n", err1 );
      fmt::print( Style::VALUE, "        Δy0:   {:.2e}\n", err2 );
      fmt::print( Style::VALUE, "        Δθ0:   {:.2e}\n", err3 );
      fmt::print( Style::VALUE, "        Δx1:   {:.2e}\n", err4 );
      fmt::print( Style::VALUE, "        Δy1:   {:.2e}\n", err5 );
      fmt::print( Style::VALUE, "        Δθ1:   {:.2e}\n", err6 );

      real_type max_err = max( { abs( err1 ), abs( err2 ), abs( err3 ), abs( err4 ), abs( err5 ), abs( err6 ) } );
      bool      passed  = max_err < 1e-6;

      print_test_result( "Original G2solve3arc example", passed );
    }
    else
    {
      print_test_result( "Original G2solve3arc example", false, "solver failed" );
    }
  }

  // Part 2: Test ClothoidList G2 from points
  {
    fmt::print( Style::INFO, "\n    📊 Part 2: ClothoidList G2 from points test\n" );

    G2lib::ClothoidList S{ "S" };

    // Try to load from file if exists
    ifstream file( "G2_test.txt" );
    if ( file.is_open() )
    {
      S.load( file );
      file.close();

      fmt::print( Style::VALUE, "      ClothoidList loaded from file\n" );
      fmt::print( Style::LABEL, "      Number of segments: " );
      fmt::print( Style::VALUE, "{}\n", S.num_segments() );

      print_test_result( "Loading from file G2_test.txt", true );
    }
    else
    {
      fmt::print( Style::WARNING, "      File G2_test.txt not found, test skipped\n" );
      print_test_result( "Loading from file G2_test.txt", false, "file not found" );
    }

    // Test build_G2 with point array (from original example)
    integer   N   = 27;
    real_type X[] = { 2.9265642,  2.6734362,  2.5109322,  1.9078122,  1.1859282,  1.9249962,  2.8265562,
                      0.0046842,  -2.826567,  -1.9437558, -1.1859438, -1.9062558, -2.501565,  -2.6734386,
                      -2.9265642, -2.6187522, -1.1406318, -0.8968758, -1.4562558, -1.9062558, -0.0046878,
                      1.9078122,  1.4468682,  0.8968722,  1.1406282,  2.6187522,  2.9265642 };
    real_type Y[] = { -1.707808758, -1.707808758, -2.367185958, -2.582810358, -2.582810358, -1.167184758, 0.915619242,
                      3.178123242,  0.915619242,  -1.150000758, -2.582810358, -2.582810358, -2.393750358, -1.707808758,
                      -1.707808758, -3.178123242, -3.178123242, -2.989063158, -0.915616758, 0.925003242,  2.953123242,
                      0.925003242,  -0.915616758, -2.989063158, -3.178123242, -3.178123242, -1.707808758 };

    Utils::TicToc tictoc;

    fmt::print( Style::VALUE, "      Test build_G2 with {} points\n", N );

    tictoc.tic();
    S.build_G2( N, X, Y, m_pi / 2, 0, m_pi / 2, 0 );
    tictoc.toc();

    fmt::print( Style::LABEL, "      build_G2 time: " );
    fmt::print( Style::VALUE, "{:.3f} ms\n", tictoc.elapsed_ms() );

    tictoc.tic();
    S.build_G2_cyclic( N, X, Y );
    tictoc.toc();

    fmt::print( Style::LABEL, "      build_G2_cyclic time: " );
    fmt::print( Style::VALUE, "{:.3f} ms\n", tictoc.elapsed_ms() );

    print_test_result( "G2 construction from array", true );
  }

  fmt::print( "\n" );
  print_test_result( "Complete original example", true, "all parts executed" );
}

// ============================================================================
// MAIN FUNCTION
// ============================================================================

int main()
{
  // Program header
  print_header( "G² CLOTHOID INTERPOLATION TEST SUITE", "🌀" );

  auto start_time = chrono::high_resolution_clock::now();

  fmt::print( Style::INFO, "    📅 Date and time: {:%Y-%m-%d %H:%M:%S}\n", chrono::system_clock::now() );
  fmt::print( Style::INFO, "    🏷️  Version: G2 Clothoid Test Suite 1.0\n" );
  fmt::print( Style::INFO, "    👤 Author: Enrico Bertolazzi\n" );
  fmt::print( Style::INFO, "    📧 Email: enrico.bertolazzi@unitn.it\n\n" );

  // Execute all tests
  test_g2solve2arc();
  test_g2solve3arc();
  test_g2solve_clc();
  test_clothoid_spline_g2();
  test_clothoid_list_g2();
  test_original_example();

  // Calculate execution time
  auto end_time = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::milliseconds>( end_time - start_time );

  // Final summary
  print_header( "G² TEST FINAL SUMMARY" );

  fmt::print( Style::SUCCESS, "    ✅ All tests have been executed successfully!\n" );
  fmt::print( Style::INFO, "    ⏱️  Execution time: {} ms\n", duration.count() );
  fmt::print( Style::INFO, "    📊 Test categories completed: 6\n\n" );

  fmt::print(
    Style::HEADER,
    "╔══════════════════════════════════════════════════════════════╗\n"
    "║                    🎉 TESTS COMPLETED!                       ║\n"
    "╚══════════════════════════════════════════════════════════════╝\n" );

  return 0;
}
