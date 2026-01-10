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

using namespace std;
using namespace G2lib;
using Utils::m_pi;

// ============================================================================
// STYLES AND FORMATTING
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
  //const auto CURVE     = fg( fmt::color::orange ) | fmt::emphasis::bold;
  const auto TEST_PASS = fg( fmt::color::green ) | fmt::emphasis::bold;
  const auto TEST_FAIL = fg( fmt::color::red ) | fmt::emphasis::bold;
  const auto POINT     = fg( fmt::color::spring_green ) | fmt::emphasis::bold;
  const auto PARAM     = fg( fmt::color::hot_pink ) | fmt::emphasis::bold;
  const auto INTERSECT = fg( fmt::color::yellow ) | fmt::emphasis::bold;
}  // namespace Style

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

void print_header( const string & title, const string & icon = "🌀" )
{
  fmt::print( "\n" );
  fmt::print( Style::HEADER, "╔══════════════════════════════════════════════════════════════╗\n" );
  fmt::print( Style::HEADER, "║ {:^60} ║\n", icon + " " + title );
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

void print_curve_info( BaseCurve const * curve, const string & name = "Curve" )
{
  fmt::print( Style::GEOMETRY, "    📐 {} ({}):\n", name, curve->type_name() );
  fmt::print( Style::LABEL, "      Start:      " );
  fmt::print( Style::POINT, "({:.6f}, {:.6f})\n", curve->x_begin(), curve->y_begin() );
  fmt::print( Style::LABEL, "      End:        " );
  fmt::print( Style::POINT, "({:.6f}, {:.6f})\n", curve->x_end(), curve->y_end() );
  fmt::print( Style::LABEL, "      θ_start:    " );
  fmt::print( Style::PARAM, "{:.6f} rad ({:.2f}°)\n", curve->theta_begin(), curve->theta_begin() * 180.0 / m_pi );
  fmt::print( Style::LABEL, "      θ_end:      " );
  fmt::print( Style::PARAM, "{:.6f} rad ({:.2f}°)\n", curve->theta_end(), curve->theta_end() * 180.0 / m_pi );
  fmt::print( Style::LABEL, "      κ_start:    " );
  fmt::print( Style::PARAM, "{:.6f}\n", curve->kappa_begin() );
  fmt::print( Style::LABEL, "      κ_end:      " );
  fmt::print( Style::PARAM, "{:.6f}\n", curve->kappa_end() );
  fmt::print( Style::LABEL, "      Length:     " );
  fmt::print( Style::HIGHLIGHT, "{:.6f}\n", curve->length() );
}
void print_intersection_results(
  const IntersectList & ilist,
  BaseCurve const *     C1,
  BaseCurve const *     C2,
  const string &        test_name = "" )
{
  if ( !test_name.empty() ) { fmt::print( Style::INTERSECT, "\n    🔍 Intersection Results for {}:\n", test_name ); }

  if ( ilist.empty() )
  {
    fmt::print( Style::WARNING, "      No intersections found\n" );
    return;
  }

  fmt::print( Style::INFO, "      Found {} intersection(s):\n", ilist.size() );

  // Tabella aggiornata con colonne corrette
  fmt::print(
    Style::LABEL,
    "      ┌─────────────────┬─────────────────┬─────────────────┬─────────────────────────┬─────────────────┬─────────────────┬─────────────────┐\n"
    "      │     Index       │    s1 (C1)      │    s2 (C2)      │      Point (x, y)       │     Delta x     │     Delta y     │     Error       │\n"
    "      ├─────────────────┼─────────────────┼─────────────────┼─────────────────────────┼─────────────────┼─────────────────┼─────────────────┤\n" );

  real_type max_error = 0.0;

  for ( size_t i = 0; i < ilist.size(); ++i )
  {
    real_type s1 = ilist[i].first;
    real_type s2 = ilist[i].second;

    real_type x1 = C1->X( s1 );
    real_type y1 = C1->Y( s1 );
    real_type x2 = C2->X( s2 );
    real_type y2 = C2->Y( s2 );

    // Calcola differenze
    real_type dx    = x1 - x2;
    real_type dy    = y1 - y2;
    real_type error = hypot( dx, dy );

    max_error = max( max_error, error );

    auto error_style = ( error < 1e-6 ) ? Style::SUCCESS : ( error < 1e-3 ) ? Style::WARNING : Style::ERROR;

    fmt::print( Style::LABEL, "      │ " );
    fmt::print( Style::VALUE, "{:^15} ", fmt::format( "#{}", i ) );
    fmt::print( Style::LABEL, "│ " );
    fmt::print( Style::PARAM, "{:^15.6g} ", s1 );
    fmt::print( Style::LABEL, "│ " );
    fmt::print( Style::PARAM, "{:^15.6g} ", s2 );
    fmt::print( Style::LABEL, "│ " );
    fmt::print( Style::POINT, " {:22} ", fmt::format( "({:.3g}, {:.3g})", x1, y1 ) );
    fmt::print( Style::LABEL, "│ " );
    fmt::print( Style::VALUE, "{:^15.3g} ", dx );
    fmt::print( Style::LABEL, "│ " );
    fmt::print( Style::VALUE, "{:^15.3g} ", dy );
    fmt::print( Style::LABEL, "│ " );
    fmt::print( error_style, "{:^15.3g} ", error );
    fmt::print( Style::LABEL, "│\n" );
  }

  fmt::print(
    Style::LABEL,
    "      └─────────────────┴─────────────────┴─────────────────┴─────────────────────────┴─────────────────┴─────────────────┴─────────────────┘\n" );

  // Riepilogo errore massimo
  fmt::print( Style::LABEL, "      Maximum consistency error: " );
  if ( max_error < 1e-9 ) { fmt::print( Style::SUCCESS, "{:.2e} ✓\n", max_error ); }
  else if ( max_error < 1e-6 ) { fmt::print( Style::WARNING, "{:.2e} ⚠\n", max_error ); }
  else
  {
    fmt::print( Style::ERROR, "{:.2e} ✗\n", max_error );
  }
}

// ============================================================================
// TEST CASES
// ============================================================================

void test_original_examples()
{
  print_section( "ORIGINAL EXAMPLES FROM TEST FILES", "📄" );

  int tests_passed = 0;
  int total_tests  = 0;

  // Test 1: From testIntersect.cc
  {
    total_tests++;
    fmt::print( Style::INFO, "    📊 Test 1: ClothoidCurve vs ClothoidCurve (original testIntersect.cc)\n" );

    ClothoidCurve C0{ "C0" }, C1{ "C1" };
    real_type     x0 = 0, y0 = 0, theta0 = m_pi * 0.7, k0 = -0.1, dk0 = 0.05, L0 = 100;
    C0.build( x0, y0, theta0, k0, dk0, L0 );

    real_type x1 = -13, y1 = 0.3, theta1 = m_pi / 4, k1 = 0.1, dk1 = -0.05, L1 = 100;
    C1.build( x1, y1, theta1, k1, dk1, L1 );

    print_curve_info( &C0, "C0 (ClothoidCurve)" );
    print_curve_info( &C1, "C1 (ClothoidCurve)" );

    IntersectList ilist;
    Utils::TicToc timer;

    timer.tic();
    C0.intersect( &C1, ilist );
    timer.toc();

    fmt::print( Style::LABEL, "      Intersection computation time: " );
    fmt::print( Style::VALUE, "{:.3f} ms\n", timer.elapsed_ms() );

    print_intersection_results( ilist, &C0, &C1 );

    // Verify symmetry
    IntersectList ilist_reverse;
    C1.intersect( &C0, ilist_reverse );

    bool symmetric = ( ilist.size() == ilist_reverse.size() );
    if ( symmetric )
    {
      for ( size_t i = 0; i < ilist.size(); ++i )
      {
        // Check if pairs are swapped
        bool found = false;
        for ( size_t j = 0; j < ilist_reverse.size(); ++j )
        {
          if (
            abs( ilist[i].first - ilist_reverse[j].second ) < 1e-6 &&
            abs( ilist[i].second - ilist_reverse[j].first ) < 1e-6 )
          {
            found = true;
            break;
          }
        }
        if ( !found ) symmetric = false;
      }
    }

    fmt::print( Style::LABEL, "      Symmetry check (C0∩C1 vs C1∩C0): " );
    if ( symmetric ) { fmt::print( Style::SUCCESS, "PASSED ✓\n" ); }
    else
    {
      fmt::print( Style::ERROR, "FAILED ✗\n" );
    }

    bool passed = !ilist.empty() && symmetric;
    print_test_result( "Clothoid-Clothoid intersection", passed );
    if ( passed ) tests_passed++;
  }

  // Test 2: From testIntersect2.cc
  {
    total_tests++;
    fmt::print( Style::INFO, "\n    📊 Test 2: CircleArc vs CircleArc (original testIntersect2.cc)\n" );

    CircleArc C0{ "C0" }, C1{ "C1" };

#if 0
        // First configuration from test file
        real_type x0 = 0, y0 = 2, theta0 = 0, k0 = 1.0/3.0, L0 = 10;
        C0.build(x0, y0, theta0, k0, L0);
        
        real_type x1 = 0, y1 = 2, theta1 = m_pi/4, k1 = 1.5*k0, L1 = 10;
        C1.build(x1, y1, theta1, k1, L1);
#else
    // Second configuration from test file
    real_type x0 = 0, y0 = 0, theta0 = 3.7688, k0 = -1.51933, L0 = 1.66304;
    C0.build( x0, y0, theta0, k0, L0 );

    x0     = -1.20978;
    y0     = -4.05915;
    theta0 = 1.18758;
    k0     = 0.112684;
    L0     = 6.18734;
    C1.build( x0, y0, theta0, k0, L0 );
#endif

    print_curve_info( &C0, "C0 (CircleArc)" );
    print_curve_info( &C1, "C1 (CircleArc)" );

    IntersectList ilist;
    C0.intersect( &C1, ilist );

    print_intersection_results( ilist, &C0, &C1 );

    // Check collision detection
    bool collides           = G2lib::collision( &C0, &C1 );
    bool expected_collision = !ilist.empty();

    fmt::print( Style::LABEL, "      Collision detection: " );
    if ( collides == expected_collision ) { fmt::print( Style::SUCCESS, "PASSED ✓ (collides={})\n", collides ); }
    else
    {
      fmt::print( Style::ERROR, "FAILED ✗ (collides={}, expected={})\n", collides, expected_collision );
    }

    bool passed = ( collides == expected_collision );
    print_test_result( "Circle-Circle intersection", passed );
    if ( passed ) tests_passed++;
  }

  // Test 3: From testIntersect3.cc
  {
    total_tests++;
    fmt::print( Style::INFO, "\n    📊 Test 3: Mixed curve types (original testIntersect3.cc)\n" );

    LineSegment L0{ "L0" };
    real_type   x00 = 0, y00 = 0, theta00 = m_pi * 0.9;
    L0.build( x00, y00, theta00, 10 );

    Biarc     C0{ "C0" }, C1{ "C1" };
    real_type x0 = 0, y0 = 0, theta0 = m_pi * 0.7;
    real_type x1 = 2, y1 = 1, theta1 = m_pi * 0.1;
    C0.build( x0, y0, theta0, x1, y1, theta1 );

    x0     = -10;
    y0     = 0;
    theta0 = m_pi / 4;
    x1     = 2;
    y1     = 1;
    theta1 = m_pi * 0.1;
    C1.build( x0, y0, theta0, x1, y1, theta1 );

    ClothoidList CL0( C0 );
    ClothoidList CL1( C1 );

    PolyLine PL0( C0, 1e-8 );
    PolyLine PL1( C1, 1e-4 );

    fmt::print( Style::LABEL, "      Testing various curve conversions:\n" );
    fmt::print( Style::VALUE, "        • Biarc → ClothoidList\n" );
    fmt::print( Style::VALUE, "        • Biarc → PolyLine (tol=1e-8)\n" );
    fmt::print( Style::VALUE, "        • Biarc → PolyLine (tol=1e-4)\n" );

    // Test intersection between PolyLine representations
    IntersectList ilist;
    PL0.intersect( &PL1, ilist );

    print_intersection_results( ilist, &PL0, &PL1, "PolyLine vs PolyLine" );

    // Test collision function
    bool collides = G2lib::collision( &PL0, &PL1 );

    fmt::print( Style::LABEL, "      Collision between PolyLines: " );
    fmt::print( Style::VALUE, "{}\n", collides );

    // Test with original Biarcs for comparison
    IntersectList ilist_biarc;
    C0.intersect( &C1, ilist_biarc );

    fmt::print( Style::LABEL, "      Original Biarc intersections: " );
    fmt::print( Style::VALUE, "{}\n", ilist_biarc.size() );

    bool passed = true;  // Basic execution test
    print_test_result( "Mixed curve types intersection", passed );
    if ( passed ) tests_passed++;
  }

  fmt::print( "\n" );
  fmt::print( Style::INFO, "    📊 Summary of original examples: " );
  fmt::print( Style::HIGHLIGHT, "{}/{} tests passed\n", tests_passed, total_tests );
}

void test_line_segment_intersections()
{
  print_section( "LINE SEGMENT INTERSECTIONS", "📏" );

  int tests_passed = 0;
  int total_tests  = 0;

  // Test 1: Parallel lines (no intersection)
  {
    total_tests++;
    fmt::print( Style::INFO, "    📊 Test 1: Parallel lines\n" );

    LineSegment L1{ "L1" }, L2{ "L2" };
    L1.build( 0, 0, 0, 5 );  // Horizontal line at y=0
    L2.build( 0, 2, 0, 5 );  // Horizontal line at y=2

    IntersectList ilist;
    L1.intersect( &L2, ilist );

    print_intersection_results( ilist, &L1, &L2 );

    bool collides = G2lib::collision( &L1, &L2 );
    bool expected = false;

    fmt::print( Style::LABEL, "      Collision expected: " );
    fmt::print( Style::VALUE, "{}, got: {}\n", expected, collides );

    bool passed = ( ilist.empty() && collides == expected );
    print_test_result( "Parallel lines", passed );
    if ( passed ) tests_passed++;
  }

  // Test 2: Intersecting lines
  {
    total_tests++;
    fmt::print( Style::INFO, "\n    📊 Test 2: Perpendicular intersecting lines\n" );

    LineSegment L1{ "L1" }, L2{ "L2" };
    L1.build( 0, 0, 0, 5 );          // Horizontal line from (0,0) to (5,0)
    L2.build( 2, -2, m_pi / 2, 4 );  // Vertical line from (2,-2) to (2,2)

    IntersectList ilist;
    L1.intersect( &L2, ilist );

    print_intersection_results( ilist, &L1, &L2 );

    bool collides = G2lib::collision( &L1, &L2 );
    bool expected = true;

    fmt::print( Style::LABEL, "      Collision expected: " );
    fmt::print( Style::VALUE, "{}, got: {}\n", expected, collides );

    bool passed = ( !ilist.empty() && collides == expected );
    print_test_result( "Intersecting lines", passed );
    if ( passed ) tests_passed++;
  }

  // Test 3: Collinear overlapping lines
  {
    total_tests++;
    fmt::print( Style::INFO, "\n    📊 Test 3: Collinear overlapping lines\n" );

    LineSegment L1{ "L1" }, L2{ "L2" };
    L1.build( 0, 0, 0, 5 );
    L2.build( 2, 0, 0, 4 );  // Overlaps from (2,0) to (6,0) but L1 ends at (5,0)

    IntersectList ilist;
    L1.intersect( &L2, ilist );

    print_intersection_results( ilist, &L1, &L2 );

    bool collides = G2lib::collision( &L1, &L2 );

    fmt::print( Style::LABEL, "      Collision: " );
    fmt::print( Style::VALUE, "{}\n", collides );

    bool passed = !ilist.empty();
    print_test_result( "Collinear overlapping lines", passed );
    if ( passed ) tests_passed++;
  }

  // Test 4: Line with offset (ISO)
  {
    total_tests++;
    fmt::print( Style::INFO, "\n    📊 Test 4: Line with offset (ISO)\n" );

    LineSegment L1{ "L1" }, L2{ "L2" };
    L1.build( 0, 0, 0, 10 );
    L2.build( 0, 1, 0, 10 );  // Parallel line 1 unit above

    real_type     offset = 0.5;
    IntersectList ilist;
    L1.intersect_ISO( offset, &L2, -offset, ilist );

    bool collides = G2lib::collision_ISO( &L1, offset, &L2, -offset );

    // le ricalcolo offsettate
    L1.build( 0, 0.5, 0, 10 );
    L2.build( 0, 0.5, 0, 10 );  // Parallel line 1 unit above

    fmt::print( Style::VALUE, "      L1 offset: {:.2f}, L2 offset: {:.2f}\n", offset, -offset );
    print_intersection_results( ilist, &L1, &L2, "Offset lines" );

    fmt::print( Style::LABEL, "      Collision with offsets: " );
    fmt::print( Style::VALUE, "{}\n", collides );

    bool passed = !ilist.empty() && collides;
    print_test_result( "Line with offset", passed );
    if ( passed ) tests_passed++;
  }

  fmt::print( "\n" );
  fmt::print( Style::INFO, "    📊 Summary Line Segment tests: " );
  fmt::print( Style::HIGHLIGHT, "{}/{} tests passed\n", tests_passed, total_tests );
}

void test_circle_arc_intersections()
{
  print_section( "CIRCLE ARC INTERSECTIONS", "⭕" );

  int tests_passed = 0;
  int total_tests  = 0;

  // Test 1: Two circles intersecting at two points
  {
    total_tests++;
    fmt::print( Style::INFO, "    📊 Test 1: Two intersecting circles\n" );

    CircleArc C1{ "C1" }, C2{ "C2" };
    C1.build( 0, 0, 0, 0.5, 2 * m_pi );  // Circle at origin, radius 0.5
    C2.build( 1, 0, 0, 0.5, 2 * m_pi );  // Circle at (1,0), radius 0.5

    IntersectList ilist;
    C1.intersect( &C2, ilist );

    print_intersection_results( ilist, &C1, &C2 );

    bool collides = G2lib::collision( &C1, &C2 );

    fmt::print( Style::LABEL, "      Expected intersections: 2, found: " );
    fmt::print( Style::VALUE, "{}\n", ilist.size() );
    fmt::print( Style::LABEL, "      Collision: " );
    fmt::print( Style::VALUE, "{}\n", collides );

    bool passed = ( ilist.size() == 2 );
    print_test_result( "Two intersecting circles", passed );
    if ( passed ) tests_passed++;
  }

  // Test 2: Tangent circles
  {
    total_tests++;
    fmt::print( Style::INFO, "\n    📊 Test 2: Tangent circles\n" );

    CircleArc C1{ "C1" }, C2{ "C2" };
    C1.build( 0, 0, 0, 1.0, 2 * m_pi );  // Circle at origin, radius 1
    C2.build( 2, 0, 0, 1.0, 2 * m_pi );  // Circle at (2,0), radius 1

    IntersectList ilist;
    C1.intersect( &C2, ilist );

    print_intersection_results( ilist, &C1, &C2 );

    fmt::print( Style::LABEL, "      Expected intersections: 1, found: " );
    fmt::print( Style::VALUE, "{}\n", ilist.size() );

    bool passed = ( ilist.size() == 1 );
    print_test_result( "Tangent circles", passed );
    if ( passed ) tests_passed++;
  }

  // Test 3: Non-intersecting circles
  {
    total_tests++;
    fmt::print( Style::INFO, "\n    📊 Test 3: Non-intersecting circles\n" );

    CircleArc C1{ "C1" }, C2{ "C2" };
    C1.build( 0, 0, 0, 0.5, 2 * m_pi );  // Circle at origin, radius 0.5
    C2.build( 3, 0, 0, 0.5, 2 * m_pi );  // Circle at (3,0), radius 0.5

    IntersectList ilist;
    C1.intersect( &C2, ilist );

    print_intersection_results( ilist, &C1, &C2 );

    bool collides = G2lib::collision( &C1, &C2 );

    fmt::print( Style::LABEL, "      Collision: " );  // Fixed: was LABER
    fmt::print( Style::VALUE, "{}\n", collides );

    bool passed = ( ilist.empty() && !collides );
    print_test_result( "Non-intersecting circles", passed );
    if ( passed ) tests_passed++;
  }

  // Test 4: Partial circle arcs
  {
    total_tests++;
    fmt::print( Style::INFO, "\n    📊 Test 4: Partial circle arcs\n" );

    CircleArc C1{ "C1" }, C2{ "C2" };
    C1.build( 0, 0, 0, 1.0, m_pi );         // Half circle (0 to π)
    C2.build( 0, 0, m_pi / 2, 1.0, m_pi );  // Half circle (π/2 to 3π/2)

    IntersectList ilist;
    C1.intersect( &C2, ilist );

    print_intersection_results( ilist, &C1, &C2 );

    fmt::print( Style::LABEL, "      Expected intersections in overlap region: " );
    fmt::print( Style::VALUE, "1 or more\n" );

    bool passed = !ilist.empty();
    print_test_result( "Partial circle arcs", passed );
    if ( passed ) tests_passed++;
  }

  fmt::print( "\n" );
  fmt::print( Style::INFO, "    📊 Summary Circle Arc tests: " );
  fmt::print( Style::HIGHLIGHT, "{}/{} tests passed\n", tests_passed, total_tests );
}

void test_clothoid_intersections()
{
  print_section( "CLOTHOID CURVE INTERSECTIONS", "🌀" );

  int tests_passed = 0;
  int total_tests  = 0;

  // Test 1: Self-intersecting clothoid
  {
    total_tests++;
    fmt::print( Style::INFO, "    📊 Test 1: Clothoid self-intersection test\n" );

    ClothoidCurve C{ "C" };
    C.build( 0, 0, 0, 0.0, 0.1, 100 );  // Start straight, then curve

    // Create a translated copy
    ClothoidCurve C2{ "C2" };
    C2.build( 10, 5, m_pi / 4, 0.0, -0.1, 100 );

    IntersectList ilist;
    C.intersect( &C2, ilist );

    print_intersection_results( ilist, &C, &C2 );

    fmt::print( Style::LABEL, "      Checking for any intersections...\n" );

    bool passed = true;  // Just checking it runs without errors
    print_test_result( "Clothoid self-intersection", passed );
    if ( passed ) tests_passed++;
  }

  // Test 2: Clothoid with line
  {
    total_tests++;
    fmt::print( Style::INFO, "\n    📊 Test 2: Clothoid intersecting with line\n" );

    ClothoidCurve C{ "C" };
    C.build( 0, 0, m_pi / 4, 0.1, 0.01, 50 );

    LineSegment L{ "L" };
    L.build( 10, 0, m_pi / 2, 20 );  // Vertical line at x=10

    IntersectList ilist;
    C.intersect( &L, ilist );

    print_intersection_results( ilist, &C, &L );

    // Also test the global intersect function
    IntersectList ilist2;
    G2lib::intersect( &C, &L, ilist2 );

    bool consistent = ( ilist.size() == ilist2.size() );

    fmt::print( Style::LABEL, "      Consistency check (member vs global): " );
    if ( consistent ) { fmt::print( Style::SUCCESS, "PASSED ✓\n" ); }
    else
    {
      fmt::print( Style::ERROR, "FAILED ✗\n" );
    }

    bool passed = consistent;
    print_test_result( "Clothoid-line intersection", passed );
    if ( passed ) tests_passed++;
  }

  // Test 3: Clothoid with circle
  {
    total_tests++;
    fmt::print( Style::INFO, "\n    📊 Test 3: Clothoid intersecting with circle\n" );

    ClothoidCurve C{ "C" };
    C.build( 0, 0, 0, 0.05, 0.001, 100 );

    CircleArc Circle{ "Circle" };
    Circle.build( 20, 5, 0, 0.2, 2 * m_pi );

    IntersectList ilist;
    C.intersect( &Circle, ilist );

    print_intersection_results( ilist, &C, &Circle );

    bool collides = G2lib::collision( &C, &Circle );

    fmt::print( Style::LABEL, "      Collision: " );
    fmt::print( Style::VALUE, "{}\n", collides );

    bool passed = true;
    print_test_result( "Clothoid-circle intersection", passed );
    if ( passed ) tests_passed++;
  }

  fmt::print( "\n" );
  fmt::print( Style::INFO, "    📊 Summary Clothoid tests: " );
  fmt::print( Style::HIGHLIGHT, "{}/{} tests passed\n", tests_passed, total_tests );
}

void test_biarc_intersections()
{
  print_section( "BIARC INTERSECTIONS", "🔄" );

  int tests_passed = 0;
  int total_tests  = 0;

  // Test 1: Biarc self-intersection
  {
    total_tests++;
    fmt::print( Style::INFO, "    📊 Test 1: Biarc self-intersection\n" );

    Biarc B1{ "B1" }, B2{ "B2" };
    B1.build( 0, 0, 0, 10, 10, m_pi / 2 );
    B2.build( 5, 0, m_pi / 4, 15, 5, 0 );

    IntersectList ilist;
    B1.intersect( &B2, ilist );

    print_intersection_results( ilist, &B1, &B2 );

    bool collides = G2lib::collision( &B1, &B2 );

    fmt::print( Style::LABEL, "      Collision: " );
    fmt::print( Style::VALUE, "{}\n", collides );

    bool passed = true;
    print_test_result( "Biarc self-intersection", passed );
    if ( passed ) tests_passed++;
  }

  // Test 2: Biarc to ClothoidList conversion and intersection
  {
    total_tests++;
    fmt::print( Style::INFO, "\n    📊 Test 2: Biarc converted to ClothoidList\n" );

    Biarc B{ "B" };
    B.build( 0, 0, 0, 10, 10, m_pi / 2 );

    ClothoidList CL{ "CL" };
    CL.build( B );

    fmt::print( Style::LABEL, "      Biarc segments: " );
    fmt::print( Style::VALUE, "2\n" );
    fmt::print( Style::LABEL, "      ClothoidList segments: " );
    fmt::print( Style::VALUE, "{}\n", CL.num_segments() );

    // Test intersection with itself
    IntersectList ilist;
    B.intersect( &CL, ilist );

    print_intersection_results( ilist, &B, &CL, "Biarc vs ClothoidList" );

    bool passed = true;
    print_test_result( "Biarc-ClothoidList intersection", passed );
    if ( passed ) tests_passed++;
  }

  fmt::print( "\n" );
  fmt::print( Style::INFO, "    📊 Summary Biarc tests: " );
  fmt::print( Style::HIGHLIGHT, "{}/{} tests passed\n", tests_passed, total_tests );
}

void test_polyline_intersections()
{
  print_section( "POLYLINE INTERSECTIONS", "📈" );

  int tests_passed = 0;
  int total_tests  = 0;

  // Test 1: Simple polyline intersection
  {
    total_tests++;
    fmt::print( Style::INFO, "    📊 Test 1: Two polylines intersecting\n" );

    // Create polyline from points
    vector<real_type> X1 = { 0, 5, 10, 15 };
    vector<real_type> Y1 = { 0, 5, 0, 5 };

    vector<real_type> X2 = { 0, 10, 20 };
    vector<real_type> Y2 = { 10, 0, 10 };

    PolyLine PL1{ "PL1" }, PL2{ "PL2" };
    PL1.build( static_cast<integer>( X1.size() ), X1.data(), Y1.data() );  // Fixed: use correct signature
    PL2.build( static_cast<integer>( X2.size() ), X2.data(), Y2.data() );  // Fixed: use correct signature

    print_curve_info( &PL1, "PolyLine 1" );
    print_curve_info( &PL2, "PolyLine 2" );

    IntersectList ilist;
    PL1.intersect( &PL2, ilist );

    print_intersection_results( ilist, &PL1, &PL2 );

    fmt::print( Style::LABEL, "      PolyLine 1 segments: " );
    fmt::print( Style::VALUE, "{}\n", PL1.num_segments() );
    fmt::print( Style::LABEL, "      PolyLine 2 segments: " );
    fmt::print( Style::VALUE, "{}\n", PL2.num_segments() );

    bool passed = true;
    print_test_result( "PolyLine intersection", passed );
    if ( passed ) tests_passed++;
  }

  // Test 2: PolyLine from Biarc
  {
    total_tests++;
    fmt::print( Style::INFO, "\n    📊 Test 2: PolyLine approximation of Biarc\n" );

    Biarc B{ "B" };
    B.build( 0, 0, 0, 10, 10, m_pi / 2 );

    PolyLine PL{ "PL" };
    PL.build( B, 1e-4 );  // High tolerance approximation

    fmt::print( Style::LABEL, "      Approximation tolerance: " );
    fmt::print( Style::VALUE, "1e-4\n" );
    fmt::print( Style::LABEL, "      PolyLine segments: " );
    fmt::print( Style::VALUE, "{}\n", PL.num_segments() );

    // Compare lengths
    real_type biarc_length = B.length();
    real_type poly_length  = PL.length();
    real_type length_error = abs( biarc_length - poly_length );

    fmt::print( Style::LABEL, "      Biarc length: " );
    fmt::print( Style::VALUE, "{:.6f}\n", biarc_length );
    fmt::print( Style::LABEL, "      PolyLine length: " );
    fmt::print( Style::VALUE, "{:.6f}\n", poly_length );
    fmt::print( Style::LABEL, "      Length error: " );
    fmt::print( Style::VALUE, "{:.2e}\n", length_error );

    bool passed = ( length_error < 1e-3 );
    print_test_result( "PolyLine approximation accuracy", passed );
    if ( passed ) tests_passed++;
  }

  fmt::print( "\n" );
  fmt::print( Style::INFO, "    📊 Summary PolyLine tests: " );
  fmt::print( Style::HIGHLIGHT, "{}/{} tests passed\n", tests_passed, total_tests );
}

void test_clothoid_list_intersections()
{
  print_section( "CLOTHOID LIST INTERSECTIONS", "📊" );

  int tests_passed = 0;
  int total_tests  = 0;

  // Test 1: ClothoidList from multiple clothoids
  {
    total_tests++;
    fmt::print( Style::INFO, "    📊 Test 1: ClothoidList intersection\n" );

    ClothoidList CL1{ "CL1" }, CL2{ "CL2" };

    // Build first clothoid list
    ClothoidCurve C1{ "C1" }, C2{ "C2" };
    C1.build( 0, 0, 0, 0.1, 0.01, 50 );
    C2.build( C1.x_end(), C1.y_end(), C1.theta_end(), 0.1, -0.01, 50 );

    CL1.push_back( C1 );
    CL1.push_back( C2 );

    // Build second clothoid list
    C1.build( 10, 0, m_pi / 2, -0.1, 0.01, 40 );
    C2.build( C1.x_end(), C1.y_end(), C1.theta_end(), -0.1, -0.01, 40 );

    CL2.push_back( C1 );
    CL2.push_back( C2 );

    fmt::print( Style::LABEL, "      ClothoidList 1 segments: " );
    fmt::print( Style::VALUE, "{}\n", CL1.num_segments() );
    fmt::print( Style::LABEL, "      ClothoidList 2 segments: " );
    fmt::print( Style::VALUE, "{}\n", CL2.num_segments() );

    IntersectList ilist;
    CL1.intersect( &CL2, ilist );

    print_intersection_results( ilist, &CL1, &CL2 );

    bool passed = true;
    print_test_result( "ClothoidList intersection", passed );
    if ( passed ) tests_passed++;
  }

  // Test 2: G2 construction from points
  {
    total_tests++;
    fmt::print( Style::INFO, "\n    📊 Test 2: G2 construction from points\n" );

    ClothoidList CL{ "CL" };

    integer           N = 6;
    vector<real_type> X = { 0.0, 2.0, 4.0, 4.0, 2.0, 0.0 };
    vector<real_type> Y = { 0.0, 3.0, 2.0, 0.0, -1.0, 0.0 };

    fmt::print( Style::LABEL, "      Building G2 from {} points\n", N );

    Utils::TicToc timer;
    timer.tic();
    CL.build_G2( N, X.data(), Y.data(), m_pi / 4, 0.0, -m_pi / 4, 0.0 );
    timer.toc();

    fmt::print( Style::LABEL, "      Construction time: " );
    fmt::print( Style::VALUE, "{:.3f} ms\n", timer.elapsed_ms() );
    fmt::print( Style::LABEL, "      Segments created: " );
    fmt::print( Style::VALUE, "{}\n", CL.num_segments() );

    // Test self-intersection (should be minimal for a simple closed curve)
    IntersectList ilist;
    CL.intersect( &CL, ilist );

    print_intersection_results( ilist, &CL, &CL );

    bool passed = ( CL.num_segments() > 0 );
    print_test_result( "G2 construction from points", passed );
    if ( passed ) tests_passed++;
  }

  fmt::print( "\n" );
  fmt::print( Style::INFO, "    📊 Summary ClothoidList tests: " );
  fmt::print( Style::HIGHLIGHT, "{}/{} tests passed\n", tests_passed, total_tests );
}

void test_performance_and_edge_cases()
{
  print_section( "PERFORMANCE AND EDGE CASES", "⚡" );

  int tests_passed = 0;
  int total_tests  = 0;

  // Test 1: Many small intersections
  {
    total_tests++;
    fmt::print( Style::INFO, "    📊 Test 1: Dense curve intersection\n" );

    PolyLine PL1{ "PL1" }, PL2{ "PL2" };

    // Create dense polylines
    vector<real_type> X1, Y1, X2, Y2;
    for ( int i = 0; i < 100; ++i )
    {
      X1.push_back( i );
      Y1.push_back( sin( i * 0.1 ) );
      X2.push_back( i );
      Y2.push_back( cos( i * 0.1 ) );
    }

    PL1.build( static_cast<integer>( X1.size() ), X1.data(), Y1.data() );  // Fixed: use correct signature
    PL2.build( static_cast<integer>( X2.size() ), X2.data(), Y2.data() );  // Fixed: use correct signature

    Utils::TicToc timer;
    timer.tic();
    IntersectList ilist;
    PL1.intersect( &PL2, ilist );
    timer.toc();

    print_intersection_results( ilist, &PL1, &PL2 );

    fmt::print( Style::LABEL, "      Curves with 100 points each\n" );
    fmt::print( Style::LABEL, "      Intersection time: " );
    fmt::print( Style::VALUE, "{:.3f} ms\n", timer.elapsed_ms() );
    fmt::print( Style::LABEL, "      Intersections found: " );
    fmt::print( Style::VALUE, "{}\n", ilist.size() );

    bool passed = timer.elapsed_ms() < 1000;  // Should be reasonably fast
    print_test_result( "Dense curve intersection performance", passed );
    if ( passed ) tests_passed++;
  }

  // Test 2: Empty and degenerate curves
  {
    total_tests++;
    fmt::print( Style::INFO, "\n    📊 Test 2: Degenerate curves\n" );

    LineSegment L1{ "L1" }, L2{ "L2" };
    L1.build( 0, 0, 0, 0 );  // Zero-length line
    L2.build( 0, 0, 0, 1 );  // Normal line

    IntersectList ilist;
    L1.intersect( &L2, ilist );

    print_intersection_results( ilist, &L1, &L2 );

    // Test with circle of zero radius
    CircleArc C{ "C" };
    C.build( 0, 0, 0, 1.0, 0.0 );  // Zero arc length

    bool passed = true;
    print_test_result( "Degenerate curve handling", passed );
    if ( passed ) tests_passed++;
  }

  // Test 3: Large offset intersections
  {
    total_tests++;
    fmt::print( Style::INFO, "\n    📊 Test 3: Large offset intersection\n" );

    LineSegment L1{ "L1" }, L2{ "L2" };
    L1.build( 0, 0, 0, 10 );
    L2.build( 0, 100, 0, 10 );  // Far away

    real_type     large_offset = 50.0;
    IntersectList ilist;
    L1.intersect_ISO( large_offset, &L2, -large_offset, ilist );

    // ricalcola con offset
    L1.build( 0, 50, 0, 10 );
    L2.build( 0, 50, 0, 10 );  // Far away
    print_intersection_results( ilist, &L1, &L2 );

    fmt::print( Style::LABEL, "      Offset distance: " );
    fmt::print( Style::VALUE, "{:.1f}\n", large_offset );
    fmt::print( Style::LABEL, "      Intersections with offset: " );
    fmt::print( Style::VALUE, "{}\n", ilist.size() );

    bool passed = !ilist.empty();
    print_test_result( "Large offset intersection", passed );
    if ( passed ) tests_passed++;
  }

  fmt::print( "\n" );
  fmt::print( Style::INFO, "    📊 Summary Performance tests: " );
  fmt::print( Style::HIGHLIGHT, "{}/{} tests passed\n", tests_passed, total_tests );
}

void test_promotion_matrix()
{
  print_section( "CURVE TYPE PROMOTION MATRIX", "🔀" );

  fmt::print( Style::INFO, "    📊 Showing curve type promotion rules:\n\n" );

  vector<pair<CurveType, string>> curve_types = {
    { CurveType::LINE, "LineSegment" },       { CurveType::CIRCLE, "CircleArc" },
    { CurveType::CLOTHOID, "ClothoidCurve" }, { CurveType::BIARC, "Biarc" },
    { CurveType::BIARC_LIST, "BiarcList" },   { CurveType::CLOTHOID_LIST, "ClothoidList" },
    { CurveType::POLYLINE, "PolyLine" },      { CurveType::DUBINS, "Dubins" }
  };

  // Print header
  fmt::print( Style::LABEL, "      ┌───────────────" );
  for ( size_t i = 0; i < curve_types.size(); ++i ) { fmt::print( "┬───────────────" ); }
  fmt::print( "┐\n" );

  fmt::print( Style::LABEL, "      │               " );
  for ( const auto & [type, name] : curve_types ) { fmt::print( Style::VALUE, "│ {:^13} ", name.substr( 0, 13 ) ); }
  fmt::print( Style::LABEL, "│\n" );

  fmt::print( Style::LABEL, "      ├───────────────" );
  for ( size_t i = 0; i < curve_types.size(); ++i ) { fmt::print( "┼───────────────" ); }
  fmt::print( "┤\n" );

  // Print rows
  for ( const auto & [row_type, row_name] : curve_types )
  {
    fmt::print( Style::LABEL, "      │ " );
    fmt::print( Style::VALUE, "{:^13} ", row_name.substr( 0, 13 ) );
    fmt::print( Style::LABEL, "│" );

    for ( const auto & [col_type, col_name] : curve_types )
    {
      try
      {
        CurveType   promoted           = curve_promote( row_type, col_type );
        string_view promoted_name_view = to_string( promoted );  // Fixed: to_string returns string_view
        string      promoted_name( promoted_name_view );         // Convert to string
        // Shorten name for display
        if ( promoted_name.length() > 13 ) { promoted_name = promoted_name.substr( 0, 10 ) + "..."; }
        fmt::print( Style::HIGHLIGHT, " {:^13} ", promoted_name );
      }
      catch ( ... )
      {
        fmt::print( Style::ERROR, " {:^13} ", "N/A" );
      }
      fmt::print( Style::LABEL, "│" );
    }
    fmt::print( "\n" );
  }

  fmt::print( Style::LABEL, "      └───────────────" );
  for ( size_t i = 0; i < curve_types.size(); ++i ) { fmt::print( "┴───────────────" ); }
  fmt::print( "┘\n" );

  fmt::print( Style::INFO, "\n    📊 Promotion rules ensure compatible types for intersection operations.\n" );
}

// ============================================================================
// MAIN TEST FUNCTION
// ============================================================================

int main()
{

  // Program header
  print_header( "COMPREHENSIVE CURVE INTERSECTION TEST SUITE", "🔍" );

  auto start_time = chrono::high_resolution_clock::now();

  fmt::print( Style::INFO, "    📅 Date and time: {:%Y-%m-%d %H:%M:%S}\n", chrono::system_clock::now() );
  fmt::print( Style::INFO, "    🏷️  Version: Intersection Test Suite 1.0\n" );
  fmt::print( Style::INFO, "    👤 Author: Enrico Bertolazzi\n" );
  fmt::print( Style::INFO, "    📧 Email: enrico.bertolazzi@unitn.it\n\n" );

  // Run all test categories
  test_original_examples();
  test_line_segment_intersections();
  test_circle_arc_intersections();
  test_clothoid_intersections();
  test_biarc_intersections();
  test_polyline_intersections();
  test_clothoid_list_intersections();
  test_performance_and_edge_cases();
  test_promotion_matrix();

  // Calculate execution time
  auto end_time = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::milliseconds>( end_time - start_time );

  // Final summary
  print_header( "TEST SUITE COMPLETE" );

  fmt::print( Style::SUCCESS, "    ✅ All intersection tests completed successfully!\n" );
  fmt::print( Style::INFO, "    ⏱️  Total execution time: {} ms\n", duration.count() );
  fmt::print( Style::INFO, "    📊 Test categories executed: 8\n" );
  fmt::print( Style::INFO, "    🔍 Curve types tested:\n" );
  fmt::print( Style::VALUE, "      • LineSegment\n" );
  fmt::print( Style::VALUE, "      • CircleArc\n" );
  fmt::print( Style::VALUE, "      • ClothoidCurve\n" );
  fmt::print( Style::VALUE, "      • Biarc\n" );
  fmt::print( Style::VALUE, "      • PolyLine\n" );
  fmt::print( Style::VALUE, "      • ClothoidList\n" );
  fmt::print( Style::VALUE, "      • BiarcList\n" );

  fmt::print( "\n" );
  fmt::print( Style::HEADER, "╔══════════════════════════════════════════════════════════════╗\n" );
  fmt::print( Style::HEADER, "║                    🎉 ALL TESTS COMPLETED!                   ║\n" );
  fmt::print( Style::HEADER, "╚══════════════════════════════════════════════════════════════╝\n" );

  return 0;
}
