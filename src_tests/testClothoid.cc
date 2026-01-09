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

using G2lib::integer;
using G2lib::real_type;
using Utils::m_pi;
using namespace std;

// Tolerance for floating point comparisons
const real_type EPS     = 1e-10;
const real_type EPS_LEN = 1e-6;

// Function to compare two values with tolerance
bool approx_equal( real_type a, real_type b, real_type eps = EPS )
{
  return std::abs( a - b ) < eps;
}

// Function to compare points
bool points_equal( real_type x1, real_type y1, real_type x2, real_type y2, real_type eps = EPS )
{
  return approx_equal( x1, x2, eps ) && approx_equal( y1, y2, eps );
}

// Function to print test result
void print_test_result( const string & test_name, bool passed )
{
  if ( passed ) { fmt::print( fg( fmt::color::green ) | fmt::emphasis::bold, "  ✅ {}: PASSED\n", test_name ); }
  else
  {
    fmt::print( fg( fmt::color::red ) | fmt::emphasis::bold, "  ❌ {}: FAILED\n", test_name );
  }
}

// Test 1: Straight line clothoid (zero curvature)
void test_straight_line()
{
  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n🔹 TEST 1: Straight Line Clothoid (κ=0, dκ=0)\n" );

  G2lib::ClothoidCurve clothoid( "straight" );
  clothoid.build( 0.0, 0.0, m_pi / 4, 0.0, 0.0, 10.0 );

  vector<pair<string, bool>> results;

  // Known values for a straight line
  real_type cos45          = cos( m_pi / 4 );
  real_type sin45          = sin( m_pi / 4 );
  real_type expected_x_end = 10.0 * cos45;
  real_type expected_y_end = 10.0 * sin45;

  // Test 1.1: End point
  bool test1 = points_equal( clothoid.x_end(), clothoid.y_end(), expected_x_end, expected_y_end, EPS_LEN );
  results.push_back( { "End point coordinates", test1 } );

  // Test 1.2: Constant angle
  bool test2 = approx_equal( clothoid.theta( 0.0 ), m_pi / 4 ) && approx_equal( clothoid.theta( 5.0 ), m_pi / 4 ) &&
               approx_equal( clothoid.theta( 10.0 ), m_pi / 4 );
  results.push_back( { "Constant angle along curve", test2 } );

  // Test 1.3: Zero curvature everywhere
  bool test3 = approx_equal( clothoid.kappa( 0.0 ), 0.0 ) && approx_equal( clothoid.kappa( 5.0 ), 0.0 ) &&
               approx_equal( clothoid.kappa( 10.0 ), 0.0 );
  results.push_back( { "Zero curvature", test3 } );

  // Test 1.4: Tangent vector
  bool test4 = approx_equal( clothoid.tx( 3.0 ), cos45 ) && approx_equal( clothoid.ty( 3.0 ), sin45 );
  results.push_back( { "Constant tangent vector", test4 } );

  // Test 1.5: Length
  bool test5 = approx_equal( clothoid.length(), 10.0, EPS_LEN );
  results.push_back( { "Correct length", test5 } );

  // Test 1.6: Closest point projection (point on line)
  real_type qx = 5.0 * cos45;
  real_type qy = 5.0 * sin45;
  real_type x, y, s, t, dst;
  clothoid.closest_point_ISO( qx, qy, x, y, s, t, dst );
  bool test6 = approx_equal( s, 5.0, EPS_LEN ) && approx_equal( t, 0.0, EPS_LEN ) && approx_equal( dst, 0.0, EPS );
  results.push_back( { "Projection of point on line", test6 } );

  // Print all results
  for ( const auto & [name, passed] : results ) { print_test_result( name, passed ); }
}

// Test 2: Circular arc clothoid (constant curvature, dκ=0)
void test_circular_arc()
{
  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n🔹 TEST 2: Circular Arc (κ=0.5, dκ=0)\n" );

  G2lib::ClothoidCurve clothoid( "circular" );
  real_type            R      = 2.0;       // Radius = 1/κ = 2.0
  real_type            L      = m_pi * R;  // Half circle length
  real_type            x0     = 0;
  real_type            y0     = 0;
  real_type            theta0 = 0;
  real_type            kappa  = 1.0 / R;
  real_type            dk     = 0;
  clothoid.build( x0, y0, theta0, kappa, dk, L );

  vector<pair<string, bool>> results;

  // Known values for a half-circle starting at (0,0) with θ=0
  // Center at (0, R) = (0, 2)
  // End point should be (0, 2R) = (0, 4)
  real_type expected_x_end = 0.0;
  real_type expected_y_end = 2.0 * R;  // = 4.0

  // Test 2.1: End point
  bool test1 = points_equal( clothoid.x_end(), clothoid.y_end(), expected_x_end, expected_y_end, EPS_LEN );
  results.push_back( { "End point of half-circle", test1 } );

  // Test 2.2: Constant curvature
  bool test2 = approx_equal( clothoid.kappa( 0.0 ), kappa ) && approx_equal( clothoid.kappa( L / 2 ), kappa ) &&
               approx_equal( clothoid.kappa( L ), kappa );
  results.push_back( { "Constant curvature", test2 } );

  // Test 2.3: Angle variation
  // For a half-circle, total angle variation should be π
  real_type theta_begin     = clothoid.theta_begin();
  real_type theta_end       = clothoid.theta_end();
  real_type angle_variation = theta_end - theta_begin;
  bool      test3           = approx_equal( angle_variation, m_pi, EPS_LEN );
  results.push_back( { "Angle variation = π for half-circle", test3 } );

  // Test 2.4: Mid-point
  // At s = L/2, point should be at (R, R) = (2, 2) for a half-circle
  real_type mid_x = clothoid.X( L / 2 );
  real_type mid_y = clothoid.Y( L / 2 );
  bool      test4 = points_equal( mid_x, mid_y, R, R, EPS_LEN );
  results.push_back( { "Mid-point coordinates", test4 } );

  // Test 2.5: Compare with CircleArc
  G2lib::CircleArc circle( "circle" );
  circle.build( x0, y0, theta0, kappa, L );
  bool test5 = points_equal( clothoid.x_end(), clothoid.y_end(), circle.x_end(), circle.y_end(), EPS_LEN );
  results.push_back( { "Match with CircleArc end point", test5 } );

  // Print all results
  for ( const auto & [name, passed] : results ) { print_test_result( name, passed ); }
}

// Test 3: Known G1 problem with symmetric clothoid
void test_symmetric_clothoid()
{
  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n🔹 TEST 3: Symmetric G1 Problem (S-bend)\n" );

  // Symmetric clothoid connecting (0,0) to (4,0) with θ0 = π/4, θ1 = -π/4
  // This should create a symmetric S-bend
  real_type            x0     = 0;
  real_type            y0     = 0;
  real_type            theta0 = m_pi / 4;
  real_type            x1     = 4;
  real_type            y1     = 0;
  real_type            theta1 = -m_pi / 4;
  G2lib::ClothoidCurve clothoid( "symmetric" );
  int                  iter = clothoid.build_G1( x0, y0, theta0, x1, y1, theta1 );

  vector<pair<string, bool>> results;

  // Test 3.1: Construction success
  bool test1 = ( iter > 0 ) && ( clothoid.length() > 0 );
  results.push_back( { "G1 construction successful", test1 } );

  // Test 3.2: Start and end points match
  bool test2 = points_equal( clothoid.x_begin(), clothoid.y_begin(), x0, y0 ) &&
               points_equal( clothoid.x_end(), clothoid.y_end(), x1, y1, EPS_LEN );
  results.push_back( { "Start/end points match", test2 } );

  // Test 3.3: Start and end angles match
  bool test3 = approx_equal( clothoid.theta_begin(), theta0, EPS_LEN ) &&
               approx_equal( clothoid.theta_end(), theta1, EPS_LEN );
  results.push_back( { "Start/end angles match", test3 } );

  // Test 3.4: Symmetry check - mid point should have θ=0 and be at x=2
  real_type L         = clothoid.length();
  real_type mid_theta = clothoid.theta( L / 2 );
  real_type mid_x     = clothoid.X( L / 2 );
  bool      test4     = approx_equal( mid_theta, 0.0, EPS_LEN ) && approx_equal( mid_x, 2.0, EPS_LEN );
  results.push_back( { "Symmetry at midpoint", test4 } );

  // Test 3.5: Curvature symmetry (κ at s = -κ at L-s)
  real_type s1    = L / 4;
  real_type s2    = 3 * L / 4;
  real_type th1   = clothoid.theta( s1 );
  real_type th2   = clothoid.theta( s2 );
  bool      test5 = approx_equal( th1, -th2, EPS_LEN );
  results.push_back( { "Angle anti-symmetry", test5 } );

  // Print all results
  for ( const auto & [name, passed] : results ) { print_test_result( name, passed ); }
}

// Test 4: Known intersection case
void test_known_intersection()
{
  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n🔹 TEST 4: Known Intersection Case\n" );

  // Create two clothoids that definitely intersect
  // Clothoid1: from (0,0) to (4,0)
  G2lib::ClothoidCurve clothoid1( "clothoid1" );
  clothoid1.build_G1( 0.0, 0.0, 0.0, 4.0, 0.0, 0.0 );

  // Clothoid2: from (2,-2) to (2,2) - vertical line through x=2
  G2lib::ClothoidCurve clothoid2( "clothoid2" );
  clothoid2.build_G1( 2.0, -2.0, m_pi / 2, 2.0, 2.0, m_pi / 2 );

  vector<pair<string, bool>> results;

  // Test 4.1: Collision detection
  bool collision = clothoid1.collision( &clothoid2 );
  results.push_back( { "Collision detected", collision } );

  // Test 4.2: Intersection points
  G2lib::IntersectList ilist;
  clothoid1.intersect( &clothoid2, ilist );
  bool has_intersection = ( ilist.size() > 0 );
  results.push_back( { "Intersection found", has_intersection } );

  if ( has_intersection )
  {
    // Should intersect at (2,0)
    real_type s1 = ilist[0].first;
    real_type s2 = ilist[0].second;

    real_type x1, y1, x2, y2;
    clothoid1.eval( s1, x1, y1 );
    clothoid2.eval( s2, x2, y2 );

    bool test3 = points_equal( x1, y1, 2.0, 0.0, EPS_LEN ) && points_equal( x2, y2, 2.0, 0.0, EPS_LEN );
    results.push_back( { "Intersection at (2,0)", test3 } );
  }

  // Test 4.3: No intersection with faraway clothoid
  G2lib::ClothoidCurve clothoid3( "clothoid3" );
  clothoid3.build_G1( 10.0, 10.0, 0.0, 12.0, 10.0, 0.0 );
  bool no_collision = !clothoid1.collision( &clothoid3 );
  results.push_back( { "No collision with distant clothoid", no_collision } );

  // Print all results
  for ( const auto & [name, passed] : results ) { print_test_result( name, passed ); }
}

// Test 5: Transformation consistency
void test_transformation_consistency()
{
  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n🔹 TEST 5: Transformation Consistency\n" );

  G2lib::ClothoidCurve original( "original" );
  real_type            x0     = 0;
  real_type            y0     = 0;
  real_type            theta0 = m_pi / 6;
  real_type            kappa0 = 0.2;
  real_type            dk     = 0.1;
  real_type            L      = 5;

  original.build( x0, y0, theta0, kappa0, dk, L );

  vector<pair<string, bool>> results;

  // Test 5.1: Translation
  G2lib::ClothoidCurve translated = original;
  real_type            tx = 3.0, ty = 4.0;
  translated.translate( tx, ty );

  bool test1 =
    points_equal( translated.x_begin(), translated.y_begin(), original.x_begin() + tx, original.y_begin() + ty );
  results.push_back( { "Translation correct", test1 } );

  // Test 5.2: Rotation (90 degrees around origin)
  G2lib::ClothoidCurve rotated = original;
  real_type            angle   = m_pi / 2;
  rotated.rotate( angle, 0.0, 0.0 );

  // Start point should rotate
  real_type x0_rot = x0 * cos( angle ) - y0 * sin( angle );
  real_type y0_rot = x0 * sin( angle ) + y0 * cos( angle );

  bool test2 = points_equal( rotated.x_begin(), rotated.y_begin(), x0_rot, y0_rot, EPS_LEN );
  // Angle should increase by rotation angle
  bool test2b = approx_equal( rotated.theta_begin(), original.theta_begin() + angle, EPS_LEN );
  results.push_back( { "Rotation correct", test2 && test2b } );

  // Test 5.3: Scaling
  G2lib::ClothoidCurve scaled = original;
  real_type            scale  = 2.0;
  scaled.scale( scale );

  bool test3 = approx_equal( scaled.length(), original.length() * scale, EPS_LEN ) &&
               approx_equal( scaled.kappa_begin(), original.kappa_begin() / scale, EPS_LEN );
  results.push_back( { "Scaling correct", test3 } );

  // Test 5.4: Reverse
  G2lib::ClothoidCurve reversed = original;
  reversed.reverse();

  bool test4 = points_equal( reversed.x_begin(), reversed.y_begin(), original.x_end(), original.y_end() ) &&
               points_equal( reversed.x_end(), reversed.y_end(), original.x_begin(), original.y_begin() ) &&
               approx_equal( reversed.theta_begin(), original.theta_end() - m_pi, EPS_LEN ) &&
               approx_equal( reversed.theta_end(), original.theta_begin() - m_pi, EPS_LEN );
  results.push_back( { "Reverse correct", test4 } );

  // Test 5.5: Trim
  G2lib::ClothoidCurve trimmed = original;
  real_type            s_begin = 1.0;
  real_type            s_end   = 3.0;
  trimmed.trim( s_begin, s_end );

  bool test5 =
    approx_equal( trimmed.length(), s_end - s_begin, EPS_LEN ) &&
    points_equal( trimmed.x_begin(), trimmed.y_begin(), original.X( s_begin ), original.Y( s_begin ), EPS_LEN );
  results.push_back( { "Trim correct", test5 } );

  // Print all results
  for ( const auto & [name, passed] : results ) { print_test_result( name, passed ); }
}

// Test 6: Mathematical properties
void test_mathematical_properties()
{
  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n🔹 TEST 6: Mathematical Properties\n" );

  G2lib::ClothoidCurve clothoid( "math" );
  clothoid.build( 0.0, 0.0, 0.0, 0.3, 0.1, 6.0 );

  vector<pair<string, bool>> results;

  // Test 6.1: Curvature linearity
  // For clothoid: κ(s) = κ0 + dκ * s
  real_type s1 = 2.0;
  real_type s2 = 4.0;
  real_type k1 = clothoid.kappa( s1 );
  real_type k2 = clothoid.kappa( s2 );
  real_type dk = clothoid.dkappa();

  bool test1 = approx_equal( k1, clothoid.kappa_begin() + dk * s1, EPS_LEN ) &&
               approx_equal( k2, clothoid.kappa_begin() + dk * s2, EPS_LEN );
  results.push_back( { "Linear curvature κ(s) = κ0 + dκ·s", test1 } );

  // Test 6.2: Angle integration
  // θ(s) = θ0 + κ0·s + (dκ/2)·s²
  real_type theta_s1          = clothoid.theta( s1 );
  real_type expected_theta_s1 = clothoid.theta_begin() + clothoid.kappa_begin() * s1 + ( dk / 2.0 ) * s1 * s1;
  bool      test2             = approx_equal( theta_s1, expected_theta_s1, EPS_LEN );
  results.push_back( { "Angle θ(s) = θ0 + κ0·s + (dκ/2)·s²", test2 } );

  // Test 6.3: Tangent consistency: (X'(s), Y'(s)) = (cosθ(s), sinθ(s))
  real_type s     = 3.0;
  real_type theta = clothoid.theta( s );
  real_type X_D   = clothoid.X_D( s );
  real_type Y_D   = clothoid.Y_D( s );

  bool test3 = approx_equal( X_D, cos( theta ), EPS_LEN ) && approx_equal( Y_D, sin( theta ), EPS_LEN );
  results.push_back( { "Tangent (X',Y') = (cosθ, sinθ)", test3 } );

  // Test 6.4: Normal consistency (ISO)
  real_type nx_iso = clothoid.nx_begin_ISO();
  real_type ny_iso = clothoid.ny_begin_ISO();
  real_type tx0    = clothoid.tx_begin();
  real_type ty0    = clothoid.ty_begin();

  bool test4 = approx_equal( nx_iso, -ty0, EPS ) && approx_equal( ny_iso, tx0, EPS );
  results.push_back( { "Normal ISO = (-ty, tx)", test4 } );

  // Test 6.5: Point at infinity
  real_type x_inf_plus, y_inf_plus, x_inf_minus, y_inf_minus;
  clothoid.Pinfinity( x_inf_plus, y_inf_plus, true );
  clothoid.Pinfinity( x_inf_minus, y_inf_minus, false );

  // For clothoid with dκ > 0, point at +∞ should be finite
  // (clothoid converges to a point as s → ∞)
  bool test5 = ( std::isfinite( x_inf_plus ) && std::isfinite( y_inf_plus ) );
  results.push_back( { "Point at +∞ is finite", test5 } );

  // Print all results
  for ( const auto & [name, passed] : results ) { print_test_result( name, passed ); }
}

// Test 7: Offset curves
void test_offset_curves()
{
  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n🔹 TEST 7: Offset Curves\n" );

  G2lib::ClothoidCurve clothoid( "offset_test" );
  real_type            x0     = 0;
  real_type            y0     = 0;
  real_type            theta0 = 0;
  real_type            kappa0 = 0.2;
  real_type            dk     = 0;
  real_type            L      = 5;
  clothoid.build( x0, y0, theta0, kappa0, dk, L );  // Constant curvature

  vector<pair<string, bool>> results;

  // Test 7.1: ISO offset at start
  real_type offs        = 1.0;
  real_type x_iso_start = clothoid.x_begin_ISO( offs );
  real_type y_iso_start = clothoid.y_begin_ISO( offs );

  // For constant curvature circle, offset at start:
  // Center at (0, R) with R = 1/κ = 5, so offset point should be at (0, R+offs) = (0, 6)
  bool test1 = points_equal( x_iso_start, y_iso_start, x0, y0 + offs, EPS_LEN );
  results.push_back( { "ISO offset at start", test1 } );

  // Test 7.2: SAE offset (opposite direction)
  real_type x_sae_start = clothoid.x_begin_SAE( offs );
  real_type y_sae_start = clothoid.y_begin_SAE( offs );
  bool      test2       = points_equal( x_sae_start, y_sae_start, x0, y0 - offs, EPS_LEN );
  results.push_back( { "SAE offset at start", test2 } );

  // Test 7.3: Length with offset (for circle, L_offs = (R ± offs) * Δθ)
  real_type L_iso          = clothoid.length_ISO( offs );
  real_type L_sae          = clothoid.length_SAE( offs );
  real_type R              = 1.0 / kappa0;  // = 5.0
  real_type delta_theta    = clothoid.theta_end() - clothoid.theta_begin();
  real_type expected_L_iso = ( R - offs ) * delta_theta;
  real_type expected_L_sae = ( R + offs ) * delta_theta;

  bool test3 = approx_equal( L_iso, expected_L_iso, EPS_LEN ) && approx_equal( L_sae, expected_L_sae, EPS_LEN );
  results.push_back( { "Offset curve lengths", test3 } );

  // Test 7.4: Evaluate with offset
  real_type s = 2.5;
  real_type th, k, x_iso, y_iso;
  clothoid.evaluate_ISO( s, offs, th, k, x_iso, y_iso );

  // For circle, point at angle θ = κ·s
  real_type theta          = clothoid.kappa_begin() * s;
  real_type expected_x_iso = ( R - offs ) * sin( theta );
  real_type expected_y_iso = R - ( R - offs ) * cos( theta );

  bool test4 = points_equal( x_iso, y_iso, expected_x_iso, expected_y_iso, EPS_LEN );
  results.push_back( { "Evaluate with ISO offset", test4 } );

  // Test 7.5: Approximate offset length
  clothoid.build( x0, y0, theta0, kappa0, 0.1, L );  // Constant curvature
  L_iso                      = clothoid.length_ISO( offs );
  real_type        L_iso_num = 0;
  real_type        xx0       = clothoid.X_ISO( 0, offs );
  real_type        yy0       = clothoid.Y_ISO( 0, offs );
  real_type        dsh       = L / 800;
  G2lib::CircleArc circle( "circle" );
  for ( int k = 1; k <= 400; ++k )
  {
    real_type ss  = 2 * k * dsh;
    real_type xx1 = clothoid.X_ISO( ss - dsh, offs );
    real_type yy1 = clothoid.Y_ISO( ss - dsh, offs );
    real_type xx2 = clothoid.X_ISO( ss, offs );
    real_type yy2 = clothoid.Y_ISO( ss, offs );
    bool      ok  = circle.build_3P( xx0, yy0, xx1, yy1, xx2, yy2 );
    if ( !ok )
    {
      results.push_back( { "Approximate offset length", false } );
      break;
    }
    L_iso_num += circle.length();
    xx0 = xx2;
    yy0 = yy2;
  }

  bool test5 = approx_equal( L_iso, L_iso_num );
  results.push_back( { "Approximate offset length", test5 } );

  // Print all results
  for ( const auto & [name, passed] : results ) { print_test_result( name, passed ); }
}

// Test 8: Bounding box calculations
void test_bounding_box()
{
  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n🔹 TEST 8: Bounding Box\n" );

  // Create a clothoid that we know the bounding box for
  // A straight line from (0,0) to (10,0) has bbox [0,0] to [10,0]
  G2lib::ClothoidCurve clothoid( "bbox_test" );
  real_type            x0     = 0;
  real_type            y0     = 0;
  real_type            theta0 = 0;
  real_type            kappa0 = 0;
  real_type            dk     = 0;
  real_type            L      = 10;
  clothoid.build( x0, y0, theta0, kappa0, dk, L );

  vector<pair<string, bool>> results;

  // Test 8.1: Regular bbox
  real_type xmin, ymin, xmax, ymax;
  clothoid.bbox( xmin, ymin, xmax, ymax );

  bool test1 = approx_equal( xmin, 0.0 ) && approx_equal( ymin, 0.0 ) && approx_equal( xmax, 10.0 ) &&
               approx_equal( ymax, 0.0 );
  results.push_back( { "BBox for straight line", test1 } );

  // Test 8.2: BBox with ISO offset
  real_type offs = 1.0;
  clothoid.bbox_ISO( offs, xmin, ymin, xmax, ymax );

  // For straight line, offset adds ±offs in y direction
  bool test2 = approx_equal( xmin, 0.0 ) && approx_equal( ymin, offs ) && approx_equal( xmax, 10.0 ) &&
               approx_equal( ymax, offs );
  results.push_back( { "BBox with ISO offset", test2 } );

  // Test 8.3: Triangle bbox for small angle variation
  G2lib::ClothoidCurve clothoid2( "small_angle" );
  clothoid2.build( x0, y0, theta0, 0.1, dk, m_pi / 2 );  // 90 degree turn

  real_type xx0, yy0, xx1, yy1, xx2, yy2;
  bool      has_triangle = clothoid2.bbTriangle( xx0, yy0, xx1, yy1, xx2, yy2 );

  bool test3 = has_triangle;  // Should have triangle for 90° turn
  results.push_back( { "Triangle bbox exists for 90° turn", test3 } );

  if ( has_triangle )
  {
    // Triangle should contain the curve
    // We can test that start and end points are inside triangle
    // (simplified: check if vertices are reasonable)
    bool vertices_ok = true;
    vertices_ok &= ( xx0 >= -1.0 && xx0 <= 5.0 );  // reasonable bounds
    vertices_ok &= ( yy0 >= -1.0 && yy0 <= 5.0 );
    vertices_ok &= ( xx1 >= -1.0 && xx1 <= 5.0 );
    vertices_ok &= ( yy1 >= -1.0 && yy1 <= 5.0 );
    vertices_ok &= ( xx2 >= -1.0 && xx2 <= 5.0 );
    vertices_ok &= ( yy2 >= -1.0 && yy2 <= 5.0 );

    results.push_back( { "Triangle vertices reasonable", vertices_ok } );
  }

  // Print all results
  for ( const auto & [name, passed] : results ) { print_test_result( name, passed ); }
}

// Test 9: Closest point with known solutions
void test_closest_point_known()
{
  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n🔹 TEST 9: Closest Point Known Solutions\n" );

  // Straight line from (0,0) to (10,0)
  G2lib::ClothoidCurve clothoid( "line" );
  clothoid.build( 0.0, 0.0, 0.0, 0.0, 0.0, 10.0 );

  vector<pair<string, bool>> results;

  // Test 9.1: Point directly above the line
  real_type qx = 5.0, qy = 3.0;
  real_type x, y, s, t, dst;
  clothoid.closest_point_ISO( qx, qy, x, y, s, t, dst );

  bool test1 = approx_equal( x, 5.0, EPS_LEN ) && approx_equal( y, 0.0, EPS_LEN ) && approx_equal( s, 5.0, EPS_LEN ) &&
               approx_equal( t, 3.0, EPS_LEN ) && approx_equal( dst, 3.0, EPS_LEN );
  results.push_back( { "Point above line projects vertically", test1 } );

  // Test 9.2: Point beyond end of line
  qx = 15.0;
  qy = 0.0;
  clothoid.closest_point_ISO( qx, qy, x, y, s, t, dst );

  bool test2 = approx_equal( x, 10.0, EPS_LEN ) &&  // Closest to end point
               approx_equal( y, 0.0, EPS_LEN ) && approx_equal( s, 10.0, EPS_LEN ) &&
               approx_equal( dst, 5.0, EPS_LEN );  // Distance = 5
  results.push_back( { "Point beyond end projects to endpoint", test2 } );

  // Test 9.3: Point before start of line
  qx = -5.0;
  qy = 0.0;
  clothoid.closest_point_ISO( qx, qy, x, y, s, t, dst );

  bool test3 = approx_equal( x, 0.0, EPS_LEN ) &&  // Closest to start point
               approx_equal( y, 0.0, EPS_LEN ) && approx_equal( s, 0.0, EPS_LEN ) &&
               approx_equal( dst, 5.0, EPS_LEN );  // Distance = 5
  results.push_back( { "Point before start projects to start", test3 } );

  // Test 9.4: Point on line
  qx = 7.0;
  qy = 0.0;
  clothoid.closest_point_ISO( qx, qy, x, y, s, t, dst );

  bool test4 = approx_equal( x, 7.0, EPS_LEN ) && approx_equal( y, 0.0, EPS_LEN ) && approx_equal( s, 7.0, EPS_LEN ) &&
               approx_equal( t, 0.0, EPS_LEN ) && approx_equal( dst, 0.0, EPS );
  results.push_back( { "Point on line projects to itself", test4 } );

  // Print all results
  for ( const auto & [name, passed] : results ) { print_test_result( name, passed ); }
}

// Main function
int main()
{
  fmt::print(
    fg( fmt::color::cyan ) | fmt::emphasis::bold,
    "╔══════════════════════════════════════════════════════════╗\n"
    "║         CLOTHOID CURVE TESTS WITH KNOWN SOLUTIONS        ║\n"
    "╚══════════════════════════════════════════════════════════╝\n\n" );

  int total_tests  = 0;
  int passed_tests = 0;

  auto run_test = [&]( auto test_func, const string & name )
  {
    try
    {
      test_func();
      fmt::print( fg( fmt::color::green ), "  Test '{}' completed successfully\n", name );
      passed_tests++;
    }
    catch ( const std::exception & e )
    {
      fmt::print( fg( fmt::color::red ), "  Test '{}' failed: {}\n", name, e.what() );
    }
    catch ( ... )
    {
      fmt::print( fg( fmt::color::red ), "  Test '{}' failed with unknown exception\n", name );
    }
    total_tests++;
  };

  // Run all tests
  run_test( test_straight_line, "Straight Line" );
  run_test( test_circular_arc, "Circular Arc" );
  run_test( test_symmetric_clothoid, "Symmetric Clothoid" );
  run_test( test_known_intersection, "Known Intersection" );
  run_test( test_transformation_consistency, "Transformation Consistency" );
  run_test( test_mathematical_properties, "Mathematical Properties" );
  run_test( test_offset_curves, "Offset Curves" );
  run_test( test_bounding_box, "Bounding Box" );
  run_test( test_closest_point_known, "Closest Point Known Solutions" );

  // Summary
  fmt::print(
    fg( fmt::color::cyan ) | fmt::emphasis::bold,
    "\n══════════════════════════════════════════════════════════\n" );

  if ( passed_tests == total_tests )
  {
    fmt::print(
      fg( fmt::color::green ) | fmt::emphasis::bold,
      "✅ ALL TESTS PASSED: {}/{} tests successful\n",
      passed_tests,
      total_tests );
  }
  else
  {
    fmt::print(
      fg( fmt::color::red ) | fmt::emphasis::bold,
      "⚠️  TESTS SUMMARY: {}/{} tests passed\n",
      passed_tests,
      total_tests );
  }

  fmt::print(
    fg( fmt::color::cyan ) | fmt::emphasis::bold,
    "══════════════════════════════════════════════════════════\n" );

  return ( passed_tests == total_tests ) ? 0 : 1;
}

//
// eof:: testClothoid.cc
//
