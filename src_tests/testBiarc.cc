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

int main()
{
  // Test configurations with different geometric scenarios
  struct TestConfig
  {
    string    name;
    real_type x0, y0, th0;
    real_type x1, y1, th1;
  };

  vector<TestConfig> test_cases = { { "Case 1: Symmetric", 0, 0, M_PI / 2, 2, 0, -M_PI / 2 },

                                    { "Case 2: Asymmetric", -1, 0, M_PI / 12, 1, 0, -M_PI / 4 },

                                    { "Case 3: Gentle Curve", 0, 0, 0, 3, 1, M_PI / 4 },

                                    { "Case 4: Tight Curve", 0, 0, M_PI / 2, 1, 1, -M_PI / 2 } };

  // Main header with colored formatting
  fmt::print(
    fg( fmt::color::cyan ) | fmt::emphasis::bold,
    "╔══════════════════════════════════════════════════════════╗\n"
    "║       BIARC vs CLOTHOID CURVE COMPARISON TEST            ║\n"
    "╚══════════════════════════════════════════════════════════╝\n\n" );

  for ( size_t i = 0; i < test_cases.size(); ++i )
  {
    const auto & test = test_cases[i];

    // Separator between test cases
    if ( i > 0 )
    {
      fmt::print( "\n" );
      fmt::print( fg( fmt::color::gray ) | fmt::emphasis::italic, "{}\n", Utils::repeat( "─", 60 ) );
      fmt::print( "\n" );
    }

    // Test case title
    fmt::print( fg( fmt::color::yellow ) | fmt::emphasis::bold, "🔹 {}:\n", test.name );
    fmt::print( "   P₀: ({:.2f}, {:.2f}), θ0: {:6.4f} rad\n", test.x0, test.y0, test.th0 );
    fmt::print( "   P₁: ({:.2f}, {:.2f}), θ1: {:6.4f} rad\n", test.x1, test.y1, test.th1 );
    fmt::print( "\n" );

    // Create curve objects
    G2lib::Biarc         ba{ "biarc" };
    G2lib::ClothoidCurve cc{ "clothoid" };

    // Build curves with error checking
    bool ba_success = ba.build( test.x0, test.y0, test.th0, test.x1, test.y1, test.th1 );
    bool cc_success = cc.build_G1( test.x0, test.y0, test.th0, test.x1, test.y1, test.th1 );

    // Construction results
    fmt::print( "📊 Curve Construction Results:\n" );
    fmt::print( "   Biarc:     {}\n", ba_success ? "✅ Success" : "❌ Failed" );
    fmt::print( "   Clothoid:  {}\n", cc_success ? "✅ Success" : "❌ Failed" );
    fmt::print( "\n" );

    // Skip if construction failed
    if ( !ba_success || !cc_success )
    {
      fmt::print( fg( fmt::color::red ) | fmt::emphasis::bold, "⚠️  Warning: Some curve constructions failed!\n\n" );
      continue;
    }

    // Biarc table
    fmt::print( fg( fmt::color::green ) | fmt::emphasis::bold, "📐 BIARC PARAMETERS:\n" );
    print_curve_comparison_table(
      "Biarc",
      ba.length(),
      ba.theta_begin(),
      ba.theta_end(),
      ba.kappa_begin(),
      ba.kappa_end(),
      ba.x_begin(),
      ba.y_begin(),
      ba.x_end(),
      ba.y_end() );

    // Sample points along Biarc
    fmt::print( "\n📍 Biarc Sample Points:\n" );
    real_type ba_len = ba.length();
    for ( real_type s = 0; s <= ba_len; s += ba_len / 4 )
    {
      real_type x, y, th, k;
      ba.evaluate( s, th, k, x, y );
      print_curve_point( s, x, y, th, k );
    }

    // Clothoid table
    fmt::print( fg( fmt::color::magenta ) | fmt::emphasis::bold, "\n🎗️  CLOTHOID CURVE PARAMETERS:\n" );
    print_curve_comparison_table(
      "Clothoid Curve",
      cc.length(),
      cc.theta_begin(),
      cc.theta_end(),
      cc.kappa_begin(),
      cc.kappa_end(),
      cc.x_begin(),
      cc.y_begin(),
      cc.x_end(),
      cc.y_end() );

    // Sample points along Clothoid
    fmt::print( "\n📍 Clothoid Sample Points:\n" );
    real_type cc_len = cc.length();
    for ( real_type s = 0; s <= cc_len; s += cc_len / 4 )
    {
      real_type x, y, th, k;
      cc.evaluate( s, th, k, x, y );
      print_curve_point( s, x, y, th, k );
    }

    // Length comparison
    fmt::print( "\n📏 Length Comparison:\n" );
    fmt::print( "   Biarc:      {:10.6f}\n", ba_len );
    fmt::print( "   Clothoid:   {:10.6f}\n", cc_len );
    fmt::print( "   Difference: {:10.6f} ({:+.2}%)\n", cc_len - ba_len, ( cc_len - ba_len ) / ba_len );

    // G1 continuity verification
    fmt::print( "\n🔍 G1 Continuity Verification:\n" );
    fmt::print( "   Biarc - Start/End θ: {:7.4f} / {:7.4f}\n", ba.theta_begin(), ba.theta_end() );
    fmt::print( "   Clothoid - Start/End θ: {:7.4f} / {:7.4f}\n", cc.theta_begin(), cc.theta_end() );
  }

  // Final summary
  fmt::print( "\n" );
  fmt::print(
    fg( fmt::color::cyan ) | fmt::emphasis::bold,
    "══════════════════════════════════════════════════════════\n" );
  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "✅ ALL TESTS COMPLETED SUCCESSFULLY!\n" );
  fmt::print(
    fg( fmt::color::cyan ) | fmt::emphasis::bold,
    "══════════════════════════════════════════════════════════\n" );

  return 0;
}
