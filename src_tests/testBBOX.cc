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
#include "Utils_string.hh"

using G2lib::real_type;
using namespace std;

// Function to print a parameter table for curves
void print_curve_comparison_table(
  const string & curve_name,
  real_type      length,
  real_type      theta_begin,
  real_type      theta_end,
  real_type      kappa_begin,
  real_type      kappa_end,
  real_type      x_begin,
  real_type      y_begin,
  real_type      x_end,
  real_type      y_end )
{
  // Table header with Unicode box drawing
  fmt::print(
    fg( fmt::color::light_blue ) | fmt::emphasis::bold,
    "\n"
    "┌───────────────────────────────────┐\n"
    "│ {:^33} │\n"
    "├─────────────────┬─────────────────┤\n",
    curve_name );

  // Curve parameters to display
  vector<pair<string, real_type>> params = { { "Length", length },   { "θ_begin", theta_begin },
                                             { "θ_end", theta_end }, { "κ_begin", kappa_begin },
                                             { "κ_end", kappa_end }, { "x_begin", x_begin },
                                             { "y_begin", y_begin }, { "x_end", x_end },
                                             { "y_end", y_end } };

  for ( const auto & [name, value] : params ) { fmt::print( "│ {:15} │ {:15.8g} │\n", name, value ); }

  fmt::print( fg( fmt::color::light_blue ) | fmt::emphasis::bold, "└─────────────────┴─────────────────┘\n" );
}

// Function to print a point on the curve with formatted output
void print_curve_point( real_type s, real_type x, real_type y, real_type theta, real_type kappa )
{
  fmt::print(
    "  s = {:6.3f} │ x = {:8.4f} │ y = {:8.4f} │ "
    "θ = {:7.4f} │ κ = {:7.4f}\n",
    s,
    x,
    y,
    theta,
    kappa );
}

// Test bounding box for a curve
void test_bbox( G2lib::BaseCurve const & curve, const string & curve_name )
{
  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n📦 BOUNDING BOX TEST for {}:\n", curve_name );

  // Test standard BBOX
  {
    real_type xmin, ymin, xmax, ymax;
    curve.bbox( xmin, ymin, xmax, ymax );
    
    fmt::print( "  Standard BBOX:\n" );
    fmt::print( "    x: [{:.8g}, {:.8g}]  width: {:.8g}\n", xmin, xmax, xmax - xmin );
    fmt::print( "    y: [{:.8g}, {:.8g}]  height: {:.8g}\n", ymin, ymax, ymax - ymin );
    
    // Check validity
    bool valid = true;
    vector<string> issues;
    
    if ( !isfinite( xmin ) || !isfinite( xmax ) || !isfinite( ymin ) || !isfinite( ymax ) ) {
      issues.push_back( "Non-finite values in BBOX" );
      valid = false;
    }
    
    if ( xmin > xmax ) {
      issues.push_back( "xmin > xmax" );
      valid = false;
    }
    
    if ( ymin > ymax ) {
      issues.push_back( "ymin > ymax" );
      valid = false;
    }
    
    if ( valid ) {
      fmt::print( fg( fmt::color::green ), "    ✅ Valid BBOX\n" );
    } else {
      fmt::print( fg( fmt::color::red ), "    ❌ Invalid BBOX:\n" );
      for ( const auto & issue : issues ) {
        fmt::print( fg( fmt::color::red ), "      - {}\n", issue );
      }
    }
  }

  // Test BBOX with ISO offset
  {
    real_type offsets[] = { 0.5, -0.5, 1.0, -1.0 };
    for ( real_type offs : offsets ) {
      real_type xmin, ymin, xmax, ymax;
      curve.bbox_ISO( offs, xmin, ymin, xmax, ymax );
      
      fmt::print( "  ISO Offset ({:+g}) BBOX:\n", offs );
      fmt::print( "    x: [{:.8g}, {:.8g}]  width: {:.8g}\n", xmin, xmax, xmax - xmin );
      fmt::print( "    y: [{:.8g}, {:.8g}]  height: {:.8g}\n", ymin, ymax, ymax - ymin );
      
      // Quick validity check
      if ( isfinite( xmin ) && isfinite( xmax ) && isfinite( ymin ) && isfinite( ymax ) &&
           xmin <= xmax && ymin <= ymax ) {
        fmt::print( fg( fmt::color::green ), "    ✅ Valid\n" );
      } else {
        fmt::print( fg( fmt::color::red ), "    ❌ Invalid\n" );
      }
    }
  }
}

// Test bounding triangle for a curve - only for ClothoidCurve
void test_bounding_triangle( G2lib::BaseCurve const & curve, const string & curve_name )
{
  fmt::print( fg( fmt::color::magenta ) | fmt::emphasis::bold, "\n🔺 BOUNDING TRIANGLE TEST for {}:\n", curve_name );

    // Test standard triangle
    {
      real_type x0, y0, x1, y1, x2, y2;
      bool has_triangle = curve.bbTriangle( x0, y0, x1, y1, x2, y2 );
      
      if ( has_triangle ) {
        fmt::print( "  Standard Triangle:\n" );
        fmt::print( "    P0: ({:.8g}, {:.8g})\n", x0, y0 );
        fmt::print( "    P1: ({:.8g}, {:.8g})\n", x1, y1 );
        fmt::print( "    P2: ({:.8g}, {:.8g})\n", x2, y2 );
        fmt::print( fg( fmt::color::green ), "    ✅ Triangle available\n" );
        
        // Test triangle with ISO offset
        for ( real_type offs : { 0.5, -0.5 } ) {
          bool has_offset_triangle = curve.bbTriangle_ISO( offs, x0, y0, x1, y1, x2, y2 );
          
          if ( has_offset_triangle ) {
            fmt::print( "  ISO Offset ({:+g}) Triangle:\n", offs );
            fmt::print( "    P0: ({:.8g}, {:.8g})\n", x0, y0 );
            fmt::print( "    P1: ({:.8g}, {:.8g})\n", x1, y1 );
            fmt::print( "    P2: ({:.8g}, {:.8g})\n", x2, y2 );
          } else {
            fmt::print( fg( fmt::color::yellow ), "  ⚠ ISO Offset ({:+g}) Triangle not available\n", offs );
          }
        }
      } else {
        fmt::print( fg( fmt::color::yellow ), "  ⚠ Triangle method not available for this curve type\n" );
      }
    }
}

// Verify that curve points are inside BBOX and triangle
void verify_containment( G2lib::BaseCurve const & curve, const string & curve_name )
{
  fmt::print( fg( fmt::color::blue ) | fmt::emphasis::bold, "\n🔍 CONTAINMENT VERIFICATION for {}:\n", curve_name );

  // Get standard BBOX
  real_type xmin, ymin, xmax, ymax;
  curve.bbox( xmin, ymin, xmax, ymax );
  
  // Get triangle if available (only for ClothoidCurve)
  real_type tx0 = 0, ty0 = 0, tx1 = 0, ty1 = 0, tx2 = 0, ty2 = 0;
  bool has_triangle = false;
  
  if ( auto clothoid = dynamic_cast<G2lib::ClothoidCurve const*>( &curve ) ) {
    has_triangle = clothoid->bbTriangle( tx0, ty0, tx1, ty1, tx2, ty2 );
  }
  
  const int n_samples = 100;
  const real_type length = curve.length();
  const real_type tolerance = 1e-12;
  
  int bbox_outside_count = 0;
  int triangle_outside_count = 0;
  vector<pair<real_type, pair<real_type, real_type>>> bbox_outside_points;
  vector<pair<real_type, pair<real_type, real_type>>> triangle_outside_points;
  
  for ( int i = 0; i <= n_samples; ++i ) {
    real_type s = ( length * i ) / n_samples;
    real_type x, y;
    curve.eval( s, x, y );
    
    // Check BBOX containment
    bool inside_bbox = ( x >= xmin - tolerance ) && ( x <= xmax + tolerance ) &&
                       ( y >= ymin - tolerance ) && ( y <= ymax + tolerance );
    
    if ( !inside_bbox ) {
      ++bbox_outside_count;
      if ( bbox_outside_count <= 3 ) { // Store first 3 outside points
        bbox_outside_points.emplace_back( s, make_pair( x, y ) );
      }
    }
    
    // Check triangle containment if available
    if ( has_triangle ) {
      // Simple barycentric coordinates test
      real_type denom = ( ( ty1 - ty2 ) * ( tx0 - tx2 ) + ( tx2 - tx1 ) * ( ty0 - ty2 ) );
      if ( abs( denom ) > 1e-12 ) {
        real_type a = ( ( ty1 - ty2 ) * ( x - tx2 ) + ( tx2 - tx1 ) * ( y - ty2 ) ) / denom;
        real_type b = ( ( ty2 - ty0 ) * ( x - tx2 ) + ( tx0 - tx2 ) * ( y - ty2 ) ) / denom;
        real_type c = 1 - a - b;
        
        bool inside_triangle = ( a >= -tolerance ) && ( b >= -tolerance ) && ( c >= -tolerance );
        if ( !inside_triangle ) {
          ++triangle_outside_count;
          if ( triangle_outside_count <= 3 ) { // Store first 3 outside points
            triangle_outside_points.emplace_back( s, make_pair( x, y ) );
          }
        }
      }
    }
  }
  
  // Report BBOX results
  if ( bbox_outside_count == 0 ) {
    fmt::print( fg( fmt::color::green ), "  ✅ All {} points are inside BBOX\n", n_samples + 1 );
  } else {
    fmt::print( fg( fmt::color::red ), "  ❌ {} of {} points are outside BBOX\n", 
                bbox_outside_count, n_samples + 1 );
    
    fmt::print( "  First few points outside BBOX:\n" );
    for ( const auto & point : bbox_outside_points ) {
      fmt::print( "    s={:.6g}: ({:.8g}, {:.8g})\n", point.first, point.second.first, point.second.second );
      fmt::print( "      BBOX: x[{:.8g}, {:.8g}], y[{:.8g}, {:.8g}]\n", xmin, xmax, ymin, ymax );
    }
  }
  
  // Report triangle results
  if ( has_triangle ) {
    if ( triangle_outside_count == 0 ) {
      fmt::print( fg( fmt::color::green ), "  ✅ All {} points are inside triangle\n", n_samples + 1 );
    } else {
      fmt::print( fg( fmt::color::red ), "  ❌ {} of {} points are outside triangle\n", 
                  triangle_outside_count, n_samples + 1 );
      
      fmt::print( "  First few points outside triangle:\n" );
      for ( const auto & point : triangle_outside_points ) {
        fmt::print( "    s={:.6g}: ({:.8g}, {:.8g})\n", point.first, point.second.first, point.second.second );
      }
    }
  } else {
    fmt::print( fg( fmt::color::yellow ), "  ⚠ Triangle containment not checked (triangle not available)\n" );
  }
}

// Test the specific problematic case from issue #55
void test_problematic_case()
{
  fmt::print( fg( fmt::color::red ) | fmt::emphasis::bold,
    "\n"
    "══════════════════════════════════════════════════════════\n"
    "🔴 TESTING PROBLEMATIC CASE FROM ISSUE #55\n"
    "══════════════════════════════════════════════════════════\n"
  );

  real_type x0 = 294.62490616135915;
  real_type y0 = 40.399218502822201;
  real_type th0 = 2.8071926535896008;
  real_type k0 = 0.40000000000205038;
  real_type dk = 0.49046153747660820;
  real_type L = 0.19160325242174281;

  fmt::print( "Parameters:\n" );
  fmt::print( "  x0:  {:.15g}\n", x0 );
  fmt::print( "  y0:  {:.15g}\n", y0 );
  fmt::print( "  th0: {:.15g}\n", th0 );
  fmt::print( "  k0:  {:.15g}\n", k0 );
  fmt::print( "  dk:  {:.15g}\n", dk );
  fmt::print( "  L:   {:.15g}\n", L );

  G2lib::ClothoidCurve clothoid( x0, y0, th0, k0, dk, L, "problematic_case" );
  
  test_bbox( clothoid, "Problematic Clothoid" );
  test_bounding_triangle( clothoid, "Problematic Clothoid" );
  verify_containment( clothoid, "Problematic Clothoid" );
}

// Test various curve types
void test_different_curve_types()
{
  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold,
    "\n══════════════════════════════════════════════════════════\n"
    "🧪 TESTING DIFFERENT CURVE TYPES\n"
    "══════════════════════════════════════════════════════════\n" );

  // Test 1: LineSegment
  {
    fmt::print( fg( fmt::color::green ), "\n📏 LineSegment Test:\n" );
    // LineSegment constructor: (x0, y0, theta0, L, name)
    G2lib::LineSegment line( 0, 0, 0, 10, "line_test" );
    test_bbox( line, "LineSegment" );
    test_bounding_triangle( line, "LineSegment" );
    verify_containment( line, "LineSegment" );
  }

  // Test 2: CircleArc
  {
    fmt::print( fg( fmt::color::green ), "\n⭕ CircleArc Test (quarter circle):\n" );
    // CircleArc constructor: (x0, y0, theta0, k, L, name)
    // For a circle of radius 1: k = 1/radius = 1, L = arc_length = radius * angle = 1 * (π/2)
    G2lib::CircleArc circle( 0, 0, 0, 1, M_PI/2, "circle_test" );
    test_bbox( circle, "CircleArc" );
    test_bounding_triangle( circle, "CircleArc" );
    verify_containment( circle, "CircleArc" );
  }

  // Test 3: ClothoidCurve
  {
    fmt::print( fg( fmt::color::green ), "\n🎗️  ClothoidCurve Test:\n" );
    // ClothoidCurve constructor: (x0, y0, th0, k0, dk, L, name)
    G2lib::ClothoidCurve clothoid( 0, 0, 0, 0.5, 0.1, 5, "clothoid_test" );
    test_bbox( clothoid, "ClothoidCurve" );
    test_bounding_triangle( clothoid, "ClothoidCurve" );
    verify_containment( clothoid, "ClothoidCurve" );
  }
}

int main()
{
  // Main header with colored formatting
  fmt::print(
    fg( fmt::color::cyan ) | fmt::emphasis::bold,
    "╔══════════════════════════════════════════════════════════╗\n"
    "║       CLOTHOID BBOX & TRIANGLE TEST SUITE                ║\n"
    "╚══════════════════════════════════════════════════════════╝\n\n"
  );

  // Test the problematic case first
  test_problematic_case();

  // Test different curve types
  test_different_curve_types();

  // Test configurations with different geometric scenarios
  struct TestConfig
  {
    string    name;
    real_type x0, y0, th0;
    real_type x1, y1, th1;
  };

  vector<TestConfig> test_cases = {
    { "Case 1: Symmetric Curve", 0, 0, M_PI / 2, 2, 0, -M_PI / 2 },
    { "Case 2: Asymmetric Curve", -1, 0, M_PI / 12, 1, 0, -M_PI / 4 },
    { "Case 3: Gentle Curve", 0, 0, 0, 3, 1, M_PI / 4 },
    { "Case 4: Tight Curve", 0, 0, M_PI / 2, 1, 1, -M_PI / 2 },
    { "Case 5: Straight Line", 0, 0, 0, 5, 0, 0 },
    { "Case 6: Quarter Circle", 0, 0, 0, 1, 1, M_PI / 2 }
  };

  for ( size_t i = 0; i < test_cases.size(); ++i )
  {
    const auto & test = test_cases[i];

    // Separator between test cases
    fmt::print( "\n" );
    fmt::print( fg( fmt::color::yellow ) | fmt::emphasis::bold,
                "══════════════════════════════════════════════════════════\n" );
    fmt::print( fg( fmt::color::yellow ) | fmt::emphasis::bold, 
                "🔹 {}:\n", test.name );
    fmt::print( fg( fmt::color::yellow ) | fmt::emphasis::bold,
                "══════════════════════════════════════════════════════════\n" );
    fmt::print( "  P₀: ({:.2f}, {:.2f}), θ0: {:6.4f} rad\n", test.x0, test.y0, test.th0 );
    fmt::print( "  P₁: ({:.2f}, {:.2f}), θ1: {:6.4f} rad\n", test.x1, test.y1, test.th1 );
    fmt::print( "\n" );

    // Create and build clothoid
    G2lib::ClothoidCurve clothoid( "test_clothoid" );
    bool cc_success = clothoid.build_G1( test.x0, test.y0, test.th0, test.x1, test.y1, test.th1 );

    if ( !cc_success ) {
      fmt::print( fg( fmt::color::red ), "❌ Failed to build clothoid\n\n" );
      continue;
    }

    // Print clothoid parameters
    fmt::print( fg( fmt::color::green ) | fmt::emphasis::bold, "📐 CLOTHOID PARAMETERS:\n" );
    print_curve_comparison_table(
      "Clothoid",
      clothoid.length(),
      clothoid.theta_begin(),
      clothoid.theta_end(),
      clothoid.kappa_begin(),
      clothoid.kappa_end(),
      clothoid.x_begin(),
      clothoid.y_begin(),
      clothoid.x_end(),
      clothoid.y_end() );

    // Sample points along Clothoid
    fmt::print( "\n📍 Clothoid Sample Points:\n" );
    real_type cc_len = clothoid.length();
    for ( real_type s = 0; s <= cc_len; s += cc_len / 4 ) {
      real_type x, y, th, k;
      clothoid.evaluate( s, th, k, x, y );
      print_curve_point( s, x, y, th, k );
    }

    // Run BBOX and triangle tests
    test_bbox( clothoid, test.name );
    test_bounding_triangle( clothoid, test.name );
    verify_containment( clothoid, test.name );
  }

  // Additional edge case tests
  fmt::print( "\n" );
  fmt::print( fg( fmt::color::magenta ) | fmt::emphasis::bold,
              "══════════════════════════════════════════════════════════\n" );
  fmt::print( fg( fmt::color::magenta ) | fmt::emphasis::bold,
              "🧪 EDGE CASE TESTS\n" );
  fmt::print( fg( fmt::color::magenta ) | fmt::emphasis::bold,
              "══════════════════════════════════════════════════════════\n" );

  // Test 1: Very short clothoid
  {
    G2lib::ClothoidCurve short_clothoid( 0, 0, 0, 100, 0, 0.001, "short_clothoid" );
    fmt::print( fg( fmt::color::cyan ), "\n📏 Very Short Clothoid (L=0.001):\n" );
    test_bbox( short_clothoid, "Short Clothoid" );
    verify_containment( short_clothoid, "Short Clothoid" );
  }

  // Test 2: High curvature clothoid
  {
    G2lib::ClothoidCurve high_curvature( 0, 0, 0, 10, 0, 1, "high_curvature" );
    fmt::print( fg( fmt::color::cyan ), "\n🌀 High Curvature Clothoid (κ=10):\n" );
    test_bbox( high_curvature, "High Curvature" );
    verify_containment( high_curvature, "High Curvature" );
  }

  // Test 3: Large offset clothoid
  {
    G2lib::ClothoidCurve base_clothoid( 0, 0, 0, 0.5, 0.1, 5, "offset_test" );
    fmt::print( fg( fmt::color::cyan ), "\n📐 Clothoid with Large Offsets:\n" );
    
    // Test with various offsets
    for ( real_type offset : { 0.0, 1.0, -1.0, 5.0, -5.0 } ) {
      real_type xmin, ymin, xmax, ymax;
      base_clothoid.bbox_ISO( offset, xmin, ymin, xmax, ymax );
      
      fmt::print( "  Offset {:+g}: x[{:.6g}, {:.6g}], y[{:.6g}, {:.6g}]\n", 
                  offset, xmin, xmax, ymin, ymax );
      
      // Verify points with offset
      const int n = 20;
      int outside = 0;
      for ( int i = 0; i <= n; ++i ) {
        real_type s = ( 5.0 * i ) / n;
        real_type x, y;
        base_clothoid.eval_ISO( s, offset, x, y );
        
        if ( !( x >= xmin - 1e-12 && x <= xmax + 1e-12 &&
                y >= ymin - 1e-12 && y <= ymax + 1e-12 ) ) {
          ++outside;
        }
      }
      
      if ( outside == 0 ) {
        fmt::print( fg( fmt::color::green ), "    ✅ All points inside\n" );
      } else {
        fmt::print( fg( fmt::color::red ), "    ❌ {} points outside\n", outside );
      }
    }
  }

  // Final summary
  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold,
    "\n"
    "══════════════════════════════════════════════════════════\n"
    "✅ ALL TESTS COMPLETED!\n"
    "══════════════════════════════════════════════════════════\n"
  );

  return 0;
}
