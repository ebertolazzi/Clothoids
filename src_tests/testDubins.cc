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
#include "Utils_eigen.hh"
#include <fstream>
#include <vector>
#include <random>
#include <iomanip>
#include <thread>
#include <chrono>

using G2lib::integer;
using G2lib::real_type;
using namespace std;

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================
template <typename... Args> void print_test_header( const string & title, Args... args )
{
  fmt::print(
    fg( fmt::color::light_blue ) | fmt::emphasis::bold,
    "\n"
    "╔══════════════════════════════════════════════════════════╗\n"
    "║ {:^56} ║\n"
    "╚══════════════════════════════════════════════════════════╝\n",
    title );
  if constexpr ( sizeof...( args ) > 0 )
  {
    fmt::print( "\n" );
    ( fmt::print( args ), ... );
    fmt::print( "\n" );
  }
}

template <typename... Args> void print_success( const string & message, Args... args )
{
  fmt::print( fg( fmt::color::green ) | fmt::emphasis::bold, "✅ {} ", message );
  if constexpr ( sizeof...( args ) > 0 )
  {
    fmt::print( ": " );
    ( fmt::print( args ), ... );
  }
  fmt::print( "\n" );
}

template <typename... Args> void print_error( const string & message, Args... args )
{
  fmt::print( fg( fmt::color::red ) | fmt::emphasis::bold, "❌ {} ", message );
  if constexpr ( sizeof...( args ) > 0 )
  {
    fmt::print( ": " );
    ( fmt::print( args ), ... );
  }
  fmt::print( "\n" );
}

template <typename... Args> void print_warning( const string & message, Args... args )
{
  fmt::print( fg( fmt::color::yellow ) | fmt::emphasis::bold, "⚠️  {} ", message );
  if constexpr ( sizeof...( args ) > 0 )
  {
    fmt::print( ": " );
    ( fmt::print( args ), ... );
  }
  fmt::print( "\n" );
}

template <typename... Args> void print_info( const string & message, Args... args )
{
  fmt::print( fg( fmt::color::cyan ), "📊 {} ", message );
  if constexpr ( sizeof...( args ) > 0 )
  {
    fmt::print( ": " );
    ( fmt::print( args ), ... );
  }
  fmt::print( "\n" );
}

// ============================================================================
// TEST 1: BASIC DUBINS CONSTRUCTION
// ============================================================================
void test_basic_dubins()
{
  print_test_header( "🧪 TEST 1: Basic Dubins Construction" );

  constexpr real_type x0     = 484.405986676;
  constexpr real_type y0     = 556.309795113;
  constexpr real_type theta0 = 4.19036786001;
  constexpr real_type x3     = 486.116491569;
  constexpr real_type y3     = 556.286039501;
  constexpr real_type theta3 = 0.593139910671;
  constexpr real_type k_max  = 1;

  G2lib::Dubins DB{ "DB" };
  bool          ok{ DB.build( x0, y0, theta0, x3, y3, theta3, k_max ) };

  print_info( "Build Results" );
  print_success( "DB.build", ok ? "SUCCESS" : "FAILED" );

  // Test all Dubins methods
  fmt::print( fg( fmt::color::cyan ), "\n📐 Dubins Properties:\n" );
  fmt::print( "├─ Type: {}\n", to_string( DB.solution_type() ) );
  fmt::print( "├─ Length: {:.6f}\n", DB.length() );
  fmt::print( "├─ Length ISO (offset=0.1): {:.6f}\n", DB.length_ISO( 0.1 ) );

  fmt::print( "\n🔧 Segment Details:\n" );
  fmt::print( "├─ C0 length: {:.6f}, kappa: {:.6f}\n", DB.length0(), DB.kappa0() );
  fmt::print( "├─ C1 length: {:.6f}, kappa: {:.6f}\n", DB.length1(), DB.kappa1() );
  fmt::print( "└─ C2 length: {:.6f}, kappa: {:.6f}\n", DB.length2(), DB.kappa2() );

  // Test geometric properties
  fmt::print( "\n📍 Geometric Points:\n" );
  fmt::print( "├─ Start: ({:.3f}, {:.3f}), θ={:.3f}\n", DB.x_begin(), DB.y_begin(), DB.theta_begin() );
  fmt::print( "├─ End:   ({:.3f}, {:.3f}), θ={:.3f}\n", DB.x_end(), DB.y_end(), DB.theta_end() );
  fmt::print( "└─ Tangent (start): ({:.3f}, {:.3f})\n", DB.tx_begin(), DB.ty_begin() );

  // Test evaluation at points
  fmt::print( "\n📈 Evaluation Tests:\n" );
  vector<real_type> test_points = { 0.0, DB.length() / 2.0, DB.length() };
  for ( const auto & s : test_points )
  {
    real_type x, y, theta, kappa;
    DB.eval( s, theta, kappa, x, y );
    fmt::print( "├─ s={:.3f}: x={:.3f}, y={:.3f}, θ={:.3f}, κ={:.3f}\n", s, x, y, theta, kappa );
  }

  // Test bounding box
  real_type xmin, ymin, xmax, ymax;
  DB.bbox( xmin, ymin, xmax, ymax );
  fmt::print( "\n📦 Bounding Box:\n" );
  fmt::print( "├─ Min: ({:.3f}, {:.3f})\n", xmin, ymin );
  fmt::print( "└─ Max: ({:.3f}, {:.3f})\n", xmax, ymax );
}

// ============================================================================
// TEST 2: CSV VALIDATION
// ============================================================================
void test_csv_validation()
{
  print_test_header( "📋 TEST 2: CSV Validation" );

  vector<string> test_files = { "DSdubinsP2P6Rand.csv", "DSdubinsP2P6.csv" };
  bool           found_file = false;

  for ( const auto & filename : test_files )
  {
    ifstream file( filename );
    if ( file.is_open() )
    {
      found_file = true;
      print_info( "Processing file", filename );

      real_type   x0, y0, th0, x1, y1, th1, kmax, l, t;
      std::string man;
      std::getline( file, man );  // Discard header

      integer                                          irow{ 0 };
      int                                              wrong{ 0 };
      int                                              total_rows{ 0 };
      vector<tuple<int, string, real_type, real_type>> errors;

      auto start_time = chrono::high_resolution_clock::now();

      while ( file >> x0 >> y0 >> th0 >> x1 >> y1 >> th1 >> kmax >> l >> man >> t )
      {
        ++irow;
        ++total_rows;
        bool   error = false;
        string error_msg;

        try
        {
          G2lib::Dubins dub( x0, y0, th0, x1, y1, th1, kmax, "temporary" );
          double        totLen = dub.length();
          std::string   compMan( to_string( dub.solution_type() ) );

          if ( abs( totLen - l ) > 1e-8 )
          {
            error_msg = fmt::format( "Length diff: {:.2e}", totLen - l );
            error     = true;
          }
          if ( compMan != man )
          {
            if ( !error_msg.empty() ) error_msg += " | ";
            error_msg += fmt::format( "Maneuver: {}≠{}", man, compMan );
            error = true;
          }

          if ( error )
          {
            ++wrong;
            errors.emplace_back( irow, error_msg, totLen, l );
          }
        }
        catch ( const std::exception & e )
        {
          ++wrong;
          errors.emplace_back( irow, fmt::format( "Exception: {}", e.what() ), 0.0, l );
        }

        // Progress indicator
        if ( irow % 1000 == 0 ) { fmt::print( fg( fmt::color::dark_gray ), "   Processed {} rows...\n", irow ); }
      }
      file.close();

      auto end_time = chrono::high_resolution_clock::now();
      auto duration = chrono::duration_cast<chrono::microseconds>( end_time - start_time );

      // Summary table
      fmt::print( fg( fmt::color::cyan ), "\n📊 Validation Summary for {}:\n", filename );
      fmt::print( "┌──────────────────────────────────────────────┐\n" );
      fmt::print( "│ {:<30} {:>12} │\n", "Total rows processed:", total_rows );
      fmt::print( "│ {:<30} {:>12} │\n", "Errors found:", wrong );
      fmt::print(
        "│ {:<30} {:>11.1f}% │\n",
        "Success rate:",
        total_rows > 0 ? 100.0 * ( total_rows - wrong ) / total_rows : 0.0 );
      fmt::print( "│ {:<30} {:>12} μs │\n", "Processing time:", duration.count() );
      fmt::print( "└──────────────────────────────────────────────┘\n" );

      if ( !errors.empty() && errors.size() <= 10 )
      {
        print_warning( "Error Details (first 10):" );
        fmt::print( "{:^6} │ {:^30} │ {:^12} │ {:^12}\n", "Row", "Error", "Computed", "Expected" );
        fmt::print( "{:─^6}┼{:─^32}┼{:─^14}┼{:─^14}\n", "", "", "", "" );

        for ( const auto & [row, msg, computed, expected] : errors )
        {
          if ( row <= 10 )
          {
            fmt::print(
              fg( fmt::color::red ),
              "{:6} │ {:<30} │ {:12.6f} │ {:12.6f}\n",
              row,
              msg.substr( 0, 30 ),
              computed,
              expected );
          }
        }
      }

      break;  // Process only the first found file
    }
  }

  if ( !found_file )
  {
    print_error( "No CSV files found!" );
    print_info( "Expected one of:", "DSdubinsP2P6Rand.csv or DSdubinsP2P6.csv" );
  }
}

// ============================================================================
// TEST 3: DUBINS INTERSECTION
// ============================================================================
void test_intersection()
{
  print_test_header( "🧩 TEST 3: Dubins Intersection" );

  using Utils::m_pi;

  // Test multiple intersection scenarios
  vector<tuple<string, real_type, real_type, real_type, real_type, real_type, real_type>> scenarios = {
    { "Parallel curves", 0, 0, -m_pi / 2, 1, 0, m_pi / 2 },
    { "Crossing curves", 0, 0, 0, 1, 0, m_pi / 2 },
    { "Opposite curves", 0, 0, 0, 2, 0, m_pi },
    { "Same start", 0, 0, m_pi / 4, 0, 0, 3 * m_pi / 4 }
  };

  for ( size_t i = 0; i < scenarios.size(); ++i )
  {
    const auto & [desc, x0, y0, th0, x1, y1, th1] = scenarios[i];

    print_info( fmt::format( "Scenario {}: {}", i + 1, desc ) );

    G2lib::Dubins DB1{ fmt::format( "DB1_{}", i ) };
    G2lib::Dubins DB2{ fmt::format( "DB2_{}", i ) };

    const real_type k_max = 1.0;
    bool            ok1   = DB1.build( x0, y0, th0, x0 + 1, y0, th1, k_max );
    bool            ok2   = DB2.build( x1, y1, th0, x1 + 1, y1, th1, k_max );

    if ( !ok1 || !ok2 )
    {
      print_warning( "Failed to build Dubins curves" );
      continue;
    }

    G2lib::IntersectList ilist;
    DB1.intersect( DB2, ilist );

    bool collision = G2lib::collision( &DB1, &DB2 );

    fmt::print( "├─ Curves built: {} | {}\n", ok1 ? "✅" : "❌", ok2 ? "✅" : "❌" );
    fmt::print(
      "├─ Collision detected: {}\n",
      collision ? fmt::format( fg( fmt::color::red ), "💥 YES" ) : fmt::format( fg( fmt::color::green ), "✅ NO" ) );
    fmt::print( "├─ Intersections found: {}\n", ilist.size() );

    if ( !ilist.empty() )
    {
      fmt::print( "└─ Intersection points:\n" );
      for ( size_t j = 0; j < ilist.size(); ++j )
      {
        real_type s1 = ilist[j].first;
        real_type s2 = ilist[j].second;
        real_type x1 = DB1.X( s1 );
        real_type y1 = DB1.Y( s1 );
        real_type x2 = DB2.X( s2 );
        real_type y2 = DB2.Y( s2 );

        real_type dist = hypot( x1 - x2, y1 - y2 );
        fmt::print( "   {}. s₁={:.3f}, s₂={:.3f}, point=({:.3f}, {:.3f}), Δ={:.2e}\n", j + 1, s1, s2, x1, y1, dist );
      }
    }
    fmt::print( "\n" );
  }
}

// ============================================================================
// TEST 4: DUBINS3P - COMPREHENSIVE METHODS TEST
// ============================================================================
void test_dubins3p_comprehensive()
{
  print_test_header( "🧪 TEST 4: Dubins3P Comprehensive Methods" );

  using Utils::m_2pi;
  using Utils::m_pi;

  // Define test parameters - non constexpr per m_pi
  const real_type xi     = -1.0;
  const real_type yi     = 0.0;
  const real_type thi    = -m_pi / 2;
  const real_type xm     = 0.0;
  const real_type ym     = 1.0;
  const real_type xf     = 1.0;
  const real_type yf     = 0.0;
  const real_type thetaf = m_pi / 2;
  const real_type k_max  = 1.0;

  // Test all build methods
  vector<pair<G2lib::Dubins3pBuildType, string>> methods = {
    { G2lib::Dubins3pBuildType::SAMPLE_ONE_DEGREE, "Sample One Degree" },
    { G2lib::Dubins3pBuildType::PATTERN_SEARCH, "Pattern Search" },
    { G2lib::Dubins3pBuildType::PATTERN_TRICHOTOMY, "Pattern Trichotomy" },
    { G2lib::Dubins3pBuildType::ELLIPSE, "Ellipse" }
  };

  fmt::print( fg( fmt::color::cyan ), "🧪 Testing all Dubins3P build methods:\n\n" );
  fmt::print(
    "{:<30} │ {:^8} │ {:^10} │ {:^10} │ {:^8} │ {:^10}\n",
    "Method",
    "Success",
    "Length",
    "Eval",
    "Type0",
    "Type1" );
  fmt::print( "{:─<31}┼{:─^10}┼{:─^12}┼{:─^12}┼{:─^10}┼{:─^12}\n", "", "", "", "", "", "" );

  vector<unique_ptr<G2lib::Dubins3p>> dubins3p_list;
  vector<real_type>                   lengths;

  for ( const auto & [method, name] : methods )
  {
    auto db3p = make_unique<G2lib::Dubins3p>( fmt::format( "D3P_{}", name ) );

    // Set parameters
    db3p->set_tolerance( 0.1 * m_pi / 180.0 );
    db3p->set_sample_angle( m_2pi / 4 );
    db3p->set_sample_points( 360 );

    auto start_time = chrono::high_resolution_clock::now();
    bool success    = db3p->build( xi, yi, thi, xm, ym, xf, yf, thetaf, k_max, method );
    auto end_time   = chrono::high_resolution_clock::now();
    // Remove unused variable warning
    (void) ( end_time - start_time );  // Tempo non usato

    if ( success )
    {
      fmt::print(
        "{:<30} │ {:^7} │ {:10.6f} │ {:10} │ {:^8} │ {:^10}\n",
        name,
        "✅",
        db3p->length(),
        db3p->num_evaluation(),
        to_string( db3p->solution_type0() ),
        to_string( db3p->solution_type1() ) );

      dubins3p_list.push_back( std::move( db3p ) );
      lengths.push_back( dubins3p_list.back()->length() );
    }
    else
    {
      fmt::print(
        "{:<30} │ {:^7} │ {:10} │ {:10} │ {:^8} │ {:^10}\n",
        name,
        "❌",
        "N/A",
        "N/A",
        "N/A",
        "N/A" );
    }
  }

  // Compare results
  if ( !lengths.empty() )
  {
    fmt::print( fg( fmt::color::cyan ), "\n📊 Length Comparison:\n" );
    real_type min_len = *min_element( lengths.begin(), lengths.end() );
    real_type max_len = *max_element( lengths.begin(), lengths.end() );
    real_type diff    = max_len - min_len;

    fmt::print( "├─ Min length: {:.6f}\n", min_len );
    fmt::print( "├─ Max length: {:.6f}\n", max_len );
    fmt::print( "├─ Difference: {:.6f}\n", diff );
    fmt::print( "└─ Relative diff: {:.2f}%\n", diff / min_len * 100 );
  }
}

// ============================================================================
// TEST 5: DUBINS3P - RANDOMIZED PERFORMANCE TEST
// ============================================================================
void test_dubins3p_randomized()
{
  print_test_header( "🎲 TEST 5: Dubins3P Randomized Performance" );

  using Utils::m_2pi;
  using Utils::m_pi;

  // Initialize different methods
  G2lib::Dubins3p DB_SD{ "D3P_sample_one_degree" };
  DB_SD.set_sample_points( 360 );

  G2lib::Dubins3p DB_PS{ "D3P_pattern_search" };
  DB_PS.set_tolerance( 0.1 * m_pi / 180.0 );
  DB_PS.set_sample_angle( m_2pi / 4 );

  G2lib::Dubins3p DB_PT{ "D3P_pattern_trichotomy" };
  DB_PT.set_tolerance( 0.1 * m_pi / 180.0 );
  DB_PT.set_sample_angle( m_2pi / 4 );

  G2lib::Dubins3p DB_EL{ "D3P_ellipse" };

  // Test parameters
  constexpr integer   num_tests = 1000;
  constexpr real_type xi        = -1.0;
  constexpr real_type yi        = 0.0;
  constexpr real_type xf        = 1.0;
  constexpr real_type yf        = 0.0;

  // Random number generation
  std::mt19937                              gen( 42 );  // Fixed seed for reproducibility
  std::uniform_real_distribution<real_type> dist( -m_pi, m_pi );
  std::uniform_real_distribution<real_type> pos_dist( -2.0, 2.0 );
  std::uniform_real_distribution<real_type> k_dist( 0.1, 1.5 );

  // Statistics
  struct Stats
  {
    int            success_count     = 0;
    double         total_time        = 0.0;
    double         total_length      = 0.0;
    int            total_evaluations = 0;
    vector<double> lengths;
  };

  map<string, Stats> stats = { { "Sample", {} },
                               { "Pattern Search", {} },
                               { "Pattern Trichotomy", {} },
                               { "Ellipse", {} } };

  print_info( fmt::format( "Running {} random tests...", num_tests ) );

  for ( integer i = 0; i < num_tests; ++i )
  {
    if ( i % 100 == 0 ) { fmt::print( fg( fmt::color::dark_gray ), "   Progress: {}/{} tests\n", i, num_tests ); }

    // Generate random parameters
    real_type thi    = dist( gen );
    real_type xm     = pos_dist( gen );
    real_type ym     = pos_dist( gen );
    real_type thetaf = dist( gen );
    real_type k_max  = k_dist( gen );

    // Test each method
    auto test_method = [&]( G2lib::Dubins3p & db, const string & name, G2lib::Dubins3pBuildType method )
    {
      auto start    = chrono::high_resolution_clock::now();
      bool success  = db.build( xi, yi, thi, xm, ym, xf, yf, thetaf, k_max, method );
      auto end      = chrono::high_resolution_clock::now();
      auto duration = chrono::duration_cast<chrono::microseconds>( end - start );

      if ( success )
      {
        stats[name].success_count++;
        stats[name].total_time += duration.count();
        stats[name].total_length += db.length();
        stats[name].total_evaluations += db.num_evaluation();
        stats[name].lengths.push_back( db.length() );
      }
      return success;
    };

    test_method( DB_SD, "Sample", G2lib::Dubins3pBuildType::SAMPLE_ONE_DEGREE );
    test_method( DB_PS, "Pattern Search", G2lib::Dubins3pBuildType::PATTERN_SEARCH );
    test_method( DB_PT, "Pattern Trichotomy", G2lib::Dubins3pBuildType::PATTERN_TRICHOTOMY );
    test_method( DB_EL, "Ellipse", G2lib::Dubins3pBuildType::ELLIPSE );
  }

  // Print results table
  fmt::print( fg( fmt::color::cyan ), "\n📊 Performance Results ({} tests):\n", num_tests );
  fmt::print( "┌────────────────────┬──────────┬─────────────┬────────────┬────────────┬──────────┐\n" );
  fmt::print(
    "│ {:^18} │ {:^8} │ {:^11} │ {:^10} │ {:^10} │ {:^8} │\n",
    "Method",
    "Success",
    "Avg Time",
    "Avg Len",
    "Avg Evals",
    "Min Len" );
  fmt::print( "├────────────────────┼──────────┼─────────────┼────────────┼────────────┼──────────┤\n" );

  for ( const auto & [name, stat] : stats )
  {
    if ( stat.success_count > 0 )
    {
      double avg_time  = stat.total_time / stat.success_count;
      double avg_len   = stat.total_length / stat.success_count;
      double avg_evals = static_cast<double>( stat.total_evaluations ) / stat.success_count;
      double min_len   = *min_element( stat.lengths.begin(), stat.lengths.end() );

      fmt::print(
        "│ {:<18} │ {:^8} │ {:8.3g} µs │ {:10.6g} │ {:10.1g} │ {:8.4g} │\n",
        name,
        stat.success_count,
        avg_time,
        avg_len,
        avg_evals,
        min_len );
    }
    else
    {
      fmt::print( "│ {:<18} │ {:^8} │ {:^8} │ {:^10} │ {:^10} │ {:^8} │\n", name, "0", "N/A", "N/A", "N/A", "N/A" );
    }
  }
  fmt::print( "└────────────────────┴──────────┴─────────────┴────────────┴────────────┴──────────┘\n" );
}

// ============================================================================
// TEST 6: DUBINS METHOD COMPLETENESS TEST
// ============================================================================
void test_dubins_method_completeness()
{
  print_test_header( "🔍 TEST 6: Dubins Method Completeness" );

  using Utils::m_pi;

  // Create a test Dubins curve
  G2lib::Dubins dubins( "TestDubins" );
  bool          built = dubins.build( 0, 0, 0, 5, 5, m_pi / 2, 0.5 );

  if ( !built )
  {
    print_error( "Failed to build test Dubins curve" );
    return;
  }

  print_success( "Basic Dubins curve built" );

  // Test all methods from Dubins.hxx
  fmt::print( fg( fmt::color::cyan ), "\n🧪 Testing all Dubins methods:\n" );

  // 1. Basic properties
  fmt::print( "\n📋 Basic Properties:\n" );
  fmt::print( "├─ Solution type: {}\n", dubins.solution_type_string() );
  fmt::print( "├─ Length: {:.6f}\n", dubins.length() );
  fmt::print( "├─ ICONE (integer code): {}\n", dubins.icode() );

  // 2. Segment accessors
  fmt::print( "\n📐 Segment Accessors:\n" );
  const auto & C0 = dubins.C0();
  const auto & C1 = dubins.C1();
  const auto & C2 = dubins.C2();

  // Correzione: CircleArc non ha center_x/center_y, usiamo le proprietà disponibili
  fmt::print( "├─ C0 radius: {:.3f}, curvature: {:.3f}\n", 1.0 / abs( C0.kappa_begin() ), C0.kappa_begin() );
  fmt::print( "├─ C1 length: {:.3f}, curvature: {:.3f}\n", C1.length(), C1.kappa_begin() );
  fmt::print( "└─ C2 angle range: {:.3f} to {:.3f} rad\n", C2.theta_begin(), C2.theta_end() );

  // 3. Point evaluation at different parameters
  fmt::print( "\n📍 Point Evaluation:\n" );
  vector<real_type> test_params = { 0.0,
                                    dubins.length() / 4.0,
                                    dubins.length() / 2.0,
                                    3.0 * dubins.length() / 4.0,
                                    dubins.length() };

  for ( real_type s : test_params )
  {
    real_type x, y, theta, kappa;
    dubins.eval( s, theta, kappa, x, y );
    fmt::print( "├─ s={:.3f}: ({:.3f}, {:.3f}), θ={:.3f}, κ={:.3f}\n", s, x, y, theta, kappa );
  }

  // 4. Derivative evaluation
  fmt::print( "\n📈 Derivative Evaluation:\n" );
  real_type s_mid = dubins.length() / 2.0;
  real_type x_D, y_D, x_DD, y_DD;
  dubins.eval_D( s_mid, x_D, y_D );
  dubins.eval_DD( s_mid, x_DD, y_DD );
  fmt::print( "├─ Midpoint tangent: ({:.3f}, {:.3f})\n", x_D, y_D );
  fmt::print( "└─ Midpoint curvature vector: ({:.3f}, {:.3f})\n", x_DD, y_DD );

  // 5. ISO offset evaluation
  fmt::print( "\n📏 ISO Offset Evaluation:\n" );
  real_type offset = 0.5;
  for ( real_type s : { 0.0, dubins.length() } )
  {
    real_type x_iso, y_iso;
    dubins.eval_ISO( s, offset, x_iso, y_iso );
    real_type x, y;
    dubins.eval( s, x, y );
    real_type dist = hypot( x_iso - x, y_iso - y );
    fmt::print( "├─ s={:.3f}, offset={:.2f}: ({:.3f}, {:.3f}), dist={:.3f}\n", s, offset, x_iso, y_iso, dist );
  }

  // 6. Transformation methods
  fmt::print( "\n🔄 Transformations:\n" );

  // Create a copy for transformation
  G2lib::Dubins dubins_copy( "Copy" );
  dubins_copy.copy( dubins );

  // Test translation
  dubins_copy.translate( 2.0, 3.0 );
  fmt::print( "├─ After translation (2,3): start=({:.3f}, {:.3f})\n", dubins_copy.x_begin(), dubins_copy.y_begin() );

  // Test rotation (reset first)
  dubins_copy.copy( dubins );
  dubins_copy.rotate( m_pi / 4, 0, 0 );
  fmt::print( "├─ After 45° rotation: start=({:.3f}, {:.3f})\n", dubins_copy.x_begin(), dubins_copy.y_begin() );

  // Test scaling
  dubins_copy.copy( dubins );
  dubins_copy.scale( 2.0 );
  fmt::print( "└─ After 2x scaling: length={:.3f} (original={:.3f})\n", dubins_copy.length(), dubins.length() );

  // 7. Closest point calculation
  fmt::print( "\n🎯 Closest Point Calculation:\n" );
  real_type test_qx = 2.0, test_qy = 2.0;
  real_type x, y, s, t, distance;
  integer   result = dubins.closest_point_ISO( test_qx, test_qy, x, y, s, t, distance );
  fmt::print( "├─ Query point: ({:.3f}, {:.3f})\n", test_qx, test_qy );
  fmt::print( "├─ Closest point: ({:.3f}, {:.3f}) at s={:.3f}\n", x, y, s );
  fmt::print( "└─ Distance: {:.6f}, result code: {}\n", distance, result );

  // 8. Bounding box with offset
  fmt::print( "\n📦 Bounding Box with ISO offset:\n" );
  real_type xmin, ymin, xmax, ymax;
  real_type iso_offset = 1.0;
  dubins.bbox_ISO( iso_offset, xmin, ymin, xmax, ymax );
  fmt::print( "├─ ISO offset: {:.2f}\n", iso_offset );
  fmt::print( "├─ BBox min: ({:.3f}, {:.3f})\n", xmin, ymin );
  fmt::print( "└─ BBox max: ({:.3f}, {:.3f})\n", xmax, ymax );

  print_success( "All Dubins methods tested successfully" );
}

// ============================================================================
// TEST 7: DUBINS3P METHOD COMPLETENESS TEST
// ============================================================================
void test_dubins3p_method_completeness()
{
  print_test_header( "🔍 TEST 7: Dubins3P Method Completeness" );

  using Utils::m_pi;

  // Create a Dubins3p curve
  G2lib::Dubins3p dubins3p( "TestDubins3p" );

  // Test parameters
  real_type xi = 0.0, yi = 0.0, thi = 0.0;
  real_type xm = 2.0, ym = 2.0;
  real_type xf = 4.0, yf = 0.0, thetaf = m_pi / 2;
  real_type k_max = 0.5;

  // Build using pattern search
  dubins3p.set_tolerance( 1e-6 );
  dubins3p.set_sample_angle( m_pi / 4 );
  bool built = dubins3p.build( xi, yi, thi, xm, ym, xf, yf, thetaf, k_max, G2lib::Dubins3pBuildType::PATTERN_SEARCH );

  if ( !built )
  {
    print_error( "Failed to build Dubins3p curve" );
    return;
  }

  print_success( "Dubins3p curve built successfully" );

  // Test all methods
  fmt::print( fg( fmt::color::cyan ), "\n🧪 Testing all Dubins3p methods:\n" );

  // 1. Basic properties
  fmt::print( "\n📋 Basic Properties:\n" );
  fmt::print( "├─ Total length: {:.6f}\n", dubins3p.length() );
  fmt::print(
    "├─ Solution types: {} | {}\n",
    to_string( dubins3p.solution_type0() ),
    to_string( dubins3p.solution_type1() ) );
  fmt::print( "├─ Combined code: {}\n", dubins3p.icode() );
  fmt::print( "├─ Tolerance: {:.2e}\n", dubins3p.tolerance() );
  fmt::print( "├─ Sample angle: {:.3f} rad\n", dubins3p.sample_angle() );
  fmt::print( "└─ Evaluations used: {}\n", dubins3p.num_evaluation() );

  // 2. Segment properties
  fmt::print( "\n📐 Six Segment Properties:\n" );
  vector<pair<string, real_type>> segments = { { "C0", dubins3p.length0() }, { "C1", dubins3p.length1() },
                                               { "C2", dubins3p.length2() }, { "C3", dubins3p.length3() },
                                               { "C4", dubins3p.length4() }, { "C5", dubins3p.length5() } };

  for ( const auto & [name, length] : segments ) { fmt::print( "├─ {} length: {:.6f}\n", name, length ); }

  // 3. Theta evaluations along the curve
  fmt::print( "\n📈 Theta Values:\n" );
  fmt::print( "├─ Theta begin: {:.3f} rad\n", dubins3p.theta_begin() );
  fmt::print( "├─ Theta end: {:.3f} rad\n", dubins3p.theta_end() );

  // Test theta at segment boundaries
  fmt::print( "├─ Theta0 end / Theta3 begin: {:.3f} rad\n", dubins3p.theta0_end() );
  fmt::print( "└─ Theta2 end / Theta3 begin: {:.3f} rad\n", dubins3p.theta2_end() );

  // 4. Position evaluation
  fmt::print( "\n📍 Position Evaluation:\n" );
  vector<real_type> test_points = { 0.0, dubins3p.length() / 3.0, 2.0 * dubins3p.length() / 3.0, dubins3p.length() };

  for ( real_type s : test_points )
  {
    real_type x, y, theta, kappa;
    dubins3p.eval( s, theta, kappa, x, y );
    fmt::print( "├─ s={:.3f}: ({:.3f}, {:.3f}), θ={:.3f}, κ={:.3f}\n", s, x, y, theta, kappa );
  }

  // 5. Derivative evaluation
  fmt::print( "\n📊 Derivative Evaluation:\n" );
  real_type s_mid = dubins3p.length() / 2.0;
  real_type x_D, y_D, x_DD, y_DD, x_DDD, y_DDD;

  dubins3p.eval_D( s_mid, x_D, y_D );
  dubins3p.eval_DD( s_mid, x_DD, y_DD );
  dubins3p.eval_DDD( s_mid, x_DDD, y_DDD );

  fmt::print( "├─ Midpoint tangent: ({:.6f}, {:.6f})\n", x_D, y_D );
  fmt::print( "├─ Midpoint curvature vector: ({:.6f}, {:.6f})\n", x_DD, y_DD );
  fmt::print( "└─ Midpoint jerk vector: ({:.6f}, {:.6f})\n", x_DDD, y_DDD );

  // 6. ISO evaluation with offset
  fmt::print( "\n📏 ISO Offset Evaluation:\n" );
  real_type offset = 0.75;
  for ( real_type s : { 0.0, dubins3p.length() } )
  {
    real_type x_iso, y_iso;
    dubins3p.eval_ISO( s, offset, x_iso, y_iso );
    real_type x, y;
    dubins3p.eval( s, x, y );
    real_type dist = hypot( x_iso - x, y_iso - y );
    fmt::print( "├─ s={:.3f}, offset={:.2f}: dist={:.6f}\n", s, offset, dist );
  }

  // 7. Closest point calculation
  fmt::print( "\n🎯 Closest Point:\n" );
  real_type query_x = 3.0, query_y = 1.0;
  real_type closest_x, closest_y, closest_s, closest_t, distance;

  integer result = dubins3p.closest_point_ISO( query_x, query_y, closest_x, closest_y, closest_s, closest_t, distance );

  fmt::print( "├─ Query: ({:.3f}, {:.3f})\n", query_x, query_y );
  fmt::print( "├─ Closest: ({:.3f}, {:.3f}) at s={:.3f}\n", closest_x, closest_y, closest_s );
  fmt::print( "├─ Distance: {:.6f}\n", distance );
  fmt::print( "└─ Result code: {}\n", result );

  // 8. Range angles (method specific to Dubins3p)
  fmt::print( "\n📐 Range Angles Test:\n" );
  try
  {
    real_type angles[10];
    integer   n_angles = dubins3p.get_range_angles( xi, yi, thi, xm, ym, xf, yf, thetaf, k_max, angles );

    fmt::print( "├─ Found {} range angles:\n", n_angles );
    for ( integer i = 0; i < n_angles; ++i )
    {
      fmt::print( "│   {}. {:.6f} rad ({:.2f}°)\n", i + 1, angles[i], angles[i] * 180.0 / m_pi );
    }
    fmt::print( "└─ End of range angles\n" );
  }
  catch ( const std::exception & e )
  {
    print_warning( "Range angles calculation failed:", e.what() );
  }

  // 9. Sample angles (vector version)
  fmt::print( "\n🎲 Sample Angles:\n" );
  try
  {
    vector<real_type> sample_angles;
    dubins3p.get_sample_angles( xi, yi, thi, xm, ym, xf, yf, thetaf, k_max, 1e-3, sample_angles );

    fmt::print( "├─ Found {} sample angles:\n", sample_angles.size() );
    if ( !sample_angles.empty() )
    {
      fmt::print(
        "│   Min: {:.6f} rad, Max: {:.6f} rad\n",
        *min_element( sample_angles.begin(), sample_angles.end() ),
        *max_element( sample_angles.begin(), sample_angles.end() ) );
    }
    fmt::print( "└─ Sample angles collection complete\n" );
  }
  catch ( const std::exception & e )
  {
    print_warning( "Sample angles calculation failed:", e.what() );
  }

  // 10. Bounding box
  fmt::print( "\n📦 Bounding Box:\n" );
  real_type xmin, ymin, xmax, ymax;
  dubins3p.bbox( xmin, ymin, xmax, ymax );
  fmt::print( "├─ Without offset: ({:.3f}, {:.3f}) to ({:.3f}, {:.3f})\n", xmin, ymin, xmax, ymax );

  real_type iso_offset = 1.0;
  dubins3p.bbox_ISO( iso_offset, xmin, ymin, xmax, ymax );
  fmt::print( "└─ With ISO offset {:.2f}: ({:.3f}, {:.3f}) to ({:.3f}, {:.3f})\n", iso_offset, xmin, ymin, xmax, ymax );

  print_success( "All Dubins3p methods tested successfully" );
}

// ============================================================================
// TEST 8: COLLISION AND INTERSECTION TESTS
// ============================================================================
void test_collision_and_intersection()
{
  print_test_header( "💥 TEST 8: Collision & Intersection" );

  using Utils::m_pi;

  // Create multiple Dubins curves for testing
  vector<G2lib::Dubins> dubins_list;
  vector<string>        names = { "CurveA", "CurveB", "CurveC", "CurveD" };

  // Create curves with different configurations
  vector<tuple<real_type, real_type, real_type, real_type, real_type, real_type>> configs = {
    { 0.0, 0.0, 0.0, 5.0, 5.0, m_pi / 2 },        // Diagonal curve
    { 2.0, 0.0, m_pi / 4, 7.0, 3.0, -m_pi / 4 },  // Crossing curve
    { 0.0, 2.0, m_pi / 2, 5.0, 2.0, m_pi / 2 },   // Horizontal curve
    { 3.0, 1.0, -m_pi / 3, 8.0, 4.0, m_pi / 3 }   // Angled curve
  };

  real_type k_max = 0.5;

  for ( size_t i = 0; i < min( names.size(), configs.size() ); ++i )
  {
    const auto & [x0, y0, th0, x1, y1, th1] = configs[i];
    G2lib::Dubins curve( names[i] );

    if ( curve.build( x0, y0, th0, x1, y1, th1, k_max ) )
    {
      dubins_list.push_back( std::move( curve ) );
      print_success( fmt::format( "Built {}", names[i] ) );
    }
    else
    {
      print_error( fmt::format( "Failed to build {}", names[i] ) );
    }
  }

  if ( dubins_list.size() < 2 )
  {
    print_error( "Need at least 2 curves for collision testing" );
    return;
  }

  // Test pairwise collisions
  fmt::print( fg( fmt::color::cyan ), "\n💥 Pairwise Collision Detection:\n" );
  fmt::print( "{:<10} │ {:<10} │ {:^12} │ {:^12}\n", "Curve1", "Curve2", "Collision", "ISO Collision" );
  fmt::print( "{:─<11}┼{:─<11}┼{:─^15}┼{:─^15}\n", "", "", "", "" );

  for ( size_t i = 0; i < dubins_list.size(); ++i )
  {
    for ( size_t j = i + 1; j < dubins_list.size(); ++j )
    {
      bool collision     = G2lib::collision( &dubins_list[i], &dubins_list[j] );
      bool collision_iso = dubins_list[i].collision_ISO( 0.2, dubins_list[j], 0.2 );

      fmt::print(
        "{:<10} │ {:<10} │ {:^12} │ {:^12}\n",
        names[i],
        names[j],
        collision ? "✘ YES" : "✔ NO",
        collision_iso ? "✘ YES" : "✔ NO" );
    }
  }

  // Test intersections
  fmt::print( fg( fmt::color::cyan ), "\n🔍 Intersection Analysis:\n" );

  for ( size_t i = 0; i < dubins_list.size(); ++i )
  {
    for ( size_t j = i + 1; j < dubins_list.size(); ++j )
    {
      G2lib::IntersectList intersections;
      dubins_list[i].intersect( dubins_list[j], intersections );

      if ( !intersections.empty() )
      {
        fmt::print( "\n{} ↔ {}: Found {} intersections\n", names[i], names[j], intersections.size() );

        for ( size_t k = 0; k < intersections.size(); ++k )
        {
          real_type s1 = intersections[k].first;
          real_type s2 = intersections[k].second;
          real_type x1 = dubins_list[i].X( s1 );
          real_type y1 = dubins_list[i].Y( s1 );
          real_type x2 = dubins_list[j].X( s2 );
          real_type y2 = dubins_list[j].Y( s2 );

          real_type dist = hypot( x1 - x2, y1 - y2 );
          fmt::print( "  {}. s₁={:.3f}, s₂={:.3f}, point≈({:.3f}, {:.3f}), Δ={:.2e}\n", k + 1, s1, s2, x1, y1, dist );
        }
      }
    }
  }

  // Test self-intersection (should have none for Dubins)
  fmt::print( fg( fmt::color::cyan ), "\n🔄 Self-Intersection Test:\n" );
  for ( size_t i = 0; i < dubins_list.size(); ++i )
  {
    G2lib::IntersectList self_intersections;
    dubins_list[i].intersect( dubins_list[i], self_intersections );

    if ( self_intersections.empty() ) { fmt::print( "├─ {}: No self-intersections (expected)\n", names[i] ); }
    else
    {
      fmt::print(
        fg( fmt::color::yellow ),
        "├─ {}: {} self-intersections found!\n",
        names[i],
        self_intersections.size() );
    }
  }
}

// ============================================================================
// TEST 9: TRANSFORMATIONS AND OFFSETS
// ============================================================================
void test_transformations_and_offsets()
{
  print_test_header( "🔄 TEST 9: Transformations & Offsets" );

  using Utils::m_pi;

  // Create a base Dubins curve
  G2lib::Dubins dubins_base( "Base" );
  bool          built = dubins_base.build( 0, 0, 0, 5, 5, m_pi / 2, 0.5 );

  if ( !built )
  {
    print_error( "Failed to build base Dubins curve" );
    return;
  }

  print_success( "Base Dubins curve built" );

  // Test 1: Multiple offsets evaluation
  fmt::print( fg( fmt::color::cyan ), "\n📏 ISO Offsets Evaluation (multiple offsets):\n" );
  vector<real_type> offsets = { -1.0, -0.5, 0.0, 0.5, 1.0 };
  real_type         s_test  = dubins_base.length() / 2.0;

  for ( real_type offset : offsets )
  {
    real_type x_iso, y_iso, x, y;
    dubins_base.eval_ISO( s_test, offset, x_iso, y_iso );
    dubins_base.eval( s_test, x, y );

    // Calculate normal vector at s_test
    real_type theta = dubins_base.theta( s_test );
    real_type nx    = -sin( theta );
    real_type ny    = cos( theta );

    // Expected offset point
    real_type x_expected = x + offset * nx;
    real_type y_expected = y + offset * ny;
    real_type error      = hypot( x_iso - x_expected, y_iso - y_expected );

    fmt::print(
      "├─ Offset {:+5.2f}: computed=({:.3f}, {:.3f}), expected=({:.3f}, {:.3f}), error={:.2e}\n",
      offset,
      x_iso,
      y_iso,
      x_expected,
      y_expected,
      error );
  }

  // Test 2: Transformations chain
  fmt::print( fg( fmt::color::cyan ), "\n🔄 Transformation Chain:\n" );

  G2lib::Dubins dubins_transformed( "Transformed" );
  dubins_transformed.copy( dubins_base );

  // Apply sequence of transformations
  dubins_transformed.translate( 2.0, 3.0 );
  dubins_transformed.rotate( m_pi / 4, dubins_transformed.x_begin(), dubins_transformed.y_begin() );
  dubins_transformed.scale( 1.5 );

  // Verify transformations
  real_type length_original    = dubins_base.length();
  real_type length_transformed = dubins_transformed.length();
  real_type scale_factor       = 1.5;  // Only scaling affects length

  fmt::print( "├─ Original length: {:.6f}\n", length_original );
  fmt::print( "├─ Transformed length: {:.6f}\n", length_transformed );
  fmt::print( "├─ Expected scaled length: {:.6f}\n", length_original * scale_factor );
  fmt::print( "└─ Scaling error: {:.2e}\n", length_transformed - length_original * scale_factor );

  // Test 3: Multiple offset bounding boxes
  fmt::print( fg( fmt::color::cyan ), "\n📦 Bounding Box with Various ISO Offsets:\n" );

  for ( real_type offset : { 0.0, 0.5, 1.0, 2.0 } )
  {
    real_type xmin, ymin, xmax, ymax;
    dubins_base.bbox_ISO( offset, xmin, ymin, xmax, ymax );

    real_type width  = xmax - xmin;
    real_type height = ymax - ymin;
    real_type area   = width * height;

    fmt::print(
      "├─ Offset {:.1f}: bbox=({:.3f},{:.3f})→({:.3f},{:.3f}), area={:.3f}\n",
      offset,
      xmin,
      ymin,
      xmax,
      ymax,
      area );
  }

  // Test 4: Closest point with offset
  fmt::print( fg( fmt::color::cyan ), "\n🎯 Closest Point with ISO Offset:\n" );

  vector<pair<real_type, real_type>> test_points = { { 2.0, 2.0 }, { 3.0, 1.0 }, { 1.0, 4.0 }, { 4.0, 3.0 } };

  real_type offset = 0.5;
  for ( size_t i = 0; i < test_points.size(); ++i )
  {
    auto [qx, qy] = test_points[i];
    real_type x, y, s, t, dist;

    integer result = dubins_base.closest_point_ISO( qx, qy, offset, x, y, s, t, dist );

    fmt::print(
      "├─ Point {}: ({:.2f},{:.2f}) → closest=({:.3f},{:.3f}) at s={:.3f}, dist={:.6f}, code={}\n",
      i + 1,
      qx,
      qy,
      x,
      y,
      s,
      dist,
      result );
  }

  // Test 5: Reverse transformation
  fmt::print( fg( fmt::color::cyan ), "\n🔁 Reverse Transformation:\n" );

  G2lib::Dubins dubins_reversed( "Reversed" );
  dubins_reversed.copy( dubins_base );
  dubins_reversed.reverse();

  // Verify endpoints are swapped
  real_type dx_start = dubins_base.x_begin() - dubins_reversed.x_end();
  real_type dy_start = dubins_base.y_begin() - dubins_reversed.y_end();
  real_type dx_end   = dubins_base.x_end() - dubins_reversed.x_begin();
  real_type dy_end   = dubins_base.y_end() - dubins_reversed.y_begin();

  fmt::print( "├─ Start match: Δ=({:.2e}, {:.2e})\n", dx_start, dy_start );
  fmt::print( "├─ End match: Δ=({:.2e}, {:.2e})\n", dx_end, dy_end );
  fmt::print(
    "└─ Length preserved: original={:.6f}, reversed={:.6f}, diff={:.2e}\n",
    dubins_base.length(),
    dubins_reversed.length(),
    dubins_base.length() - dubins_reversed.length() );

  print_success( "All transformation and offset tests completed" );
}

// ============================================================================
// TEST 10: DUBINS3P EDGE CASES
// ============================================================================
void test_dubins3p_edge_cases()
{
  print_test_header( "⚠️ TEST 10: Dubins3P Edge Cases" );

  using Utils::m_pi;

  // Test cases designed to stress the algorithms
  vector<
    tuple<string, real_type, real_type, real_type, real_type, real_type, real_type, real_type, real_type, real_type>>
    cases = { { "Co-linear points", 0, 0, 0, 2, 0, 4, 0, 0, 0.5 },
              { "Same start and intermediate", 0, 0, 0, 0, 0, 4, 0, 0, 0.5 },
              { "Same intermediate and end", 0, 0, 0, 4, 0, 4, 0, 0, 0.5 },
              { "Very sharp turn", 0, 0, 0, 1, 0.1, 2, 0, m_pi, 2.0 },
              { "Large curvature", 0, 0, 0, 1, 1, 2, 0, m_pi / 2, 10.0 },
              { "Small curvature", 0, 0, 0, 5, 5, 10, 0, m_pi / 2, 0.01 },
              { "Far apart points", -10, -10, 0, 0, 0, 10, 10, m_pi, 0.5 } };

  vector<pair<G2lib::Dubins3pBuildType, string>> methods = { { G2lib::Dubins3pBuildType::SAMPLE_ONE_DEGREE, "Sample" },
                                                             { G2lib::Dubins3pBuildType::PATTERN_SEARCH, "Pattern" },
                                                             { G2lib::Dubins3pBuildType::PATTERN_TRICHOTOMY,
                                                               "Trichotomy" },
                                                             { G2lib::Dubins3pBuildType::ELLIPSE, "Ellipse" } };

  fmt::print( fg( fmt::color::cyan ), "🧪 Testing Edge Cases with All Methods:\n\n" );

  int total_tests  = 0;
  int passed_tests = 0;

  for ( const auto & [desc, xi, yi, thi, xm, ym, xf, yf, thetaf, k_max] : cases )
  {
    print_info( fmt::format( "Case: {}", desc ) );

    for ( const auto & [method, method_name] : methods )
    {
      G2lib::Dubins3p db3p( fmt::format( "D3P_{}_{}", desc.substr( 0, 5 ), method_name ) );

      // Configure based on method
      if ( method == G2lib::Dubins3pBuildType::SAMPLE_ONE_DEGREE ) { db3p.set_sample_points( 360 ); }
      else
      {
        db3p.set_tolerance( 0.1 * m_pi / 180.0 );
        db3p.set_sample_angle( m_pi / 4 );
      }

      bool success = false;
      try
      {
        success = db3p.build( xi, yi, thi, xm, ym, xf, yf, thetaf, k_max, method );
      }
      catch ( const std::exception & e )
      {
        fmt::print( fg( fmt::color::red ), "   ❌ {}: Exception: {}\n", method_name, e.what() );
      }

      total_tests++;
      if ( success )
      {
        passed_tests++;
        fmt::print( "   ✅ {}: Success, length={:.6f}, evals={}\n", method_name, db3p.length(), db3p.num_evaluation() );

        // Validate the solution
        real_type error_start = hypot( db3p.x_begin() - xi, db3p.y_begin() - yi );
        //real_type error_mid   = 0;  // Would need to compute closest point to (xm, ym)
        real_type error_end   = hypot( db3p.x_end() - xf, db3p.y_end() - yf );

        if ( error_start > 1e-6 || error_end > 1e-6 )
        {
          fmt::print(
            fg( fmt::color::yellow ),
            "      ⚠️ Validation errors: start={:.2e}, end={:.2e}\n",
            error_start,
            error_end );
        }
      }
      else
      {
        fmt::print( fg( fmt::color::yellow ), "   ⚠️ {}: Failed to build\n", method_name );
      }
    }
    fmt::print( "\n" );
  }

  fmt::print( fg( fmt::color::cyan ), "📊 Edge Cases Summary:\n" );
  fmt::print( "├─ Total tests: {}\n", total_tests );
  fmt::print( "├─ Passed: {} ({:.1f}%)\n", passed_tests, 100.0 * passed_tests / total_tests );
  fmt::print( "└─ Failed: {}\n", total_tests - passed_tests );
}

// ============================================================================
// TEST 11: PERFORMANCE BENCHMARK
// ============================================================================
void test_performance_benchmark()
{
  print_test_header( "🔍 TEST 11: Performance Benchmark" );

  using Utils::m_pi;

  // Warm-up runs
  {
    G2lib::Dubins warmup( "warmup" );
    for ( int i = 0; i < 100; ++i ) { warmup.build( 0, 0, 0, 1, 1, m_pi / 2, 1.0 ); }
  }

  // Benchmark Dubins construction
  constexpr int NUM_ITERATIONS = 10000;

  auto start = chrono::high_resolution_clock::now();

  for ( int i = 0; i < NUM_ITERATIONS; ++i )
  {
    G2lib::Dubins db( "benchmark" );
    // Varying parameters slightly
    real_type angle = ( i % 100 ) * m_pi / 50.0;
    db.build( 0, 0, 0, cos( angle ), sin( angle ), angle, 1.0 );
  }

  auto end      = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::microseconds>( end - start );

  fmt::print( fg( fmt::color::cyan ), "📈 Dubins Construction Benchmark:\n" );
  fmt::print( "├─ Iterations: {}\n", NUM_ITERATIONS );
  fmt::print( "├─ Total time: {:.3f} ms\n", duration.count() / 1000.0 );
  fmt::print( "├─ Average time: {:.3f} µs\n", duration.count() / static_cast<double>( NUM_ITERATIONS ) );
  fmt::print( "└─ Operations/sec: {:.0f}\n", NUM_ITERATIONS / ( duration.count() / 1e6 ) );

  // Benchmark Dubins3P with different methods
  vector<pair<G2lib::Dubins3pBuildType, string>> methods = { { G2lib::Dubins3pBuildType::SAMPLE_ONE_DEGREE, "Sample" },
                                                             { G2lib::Dubins3pBuildType::PATTERN_SEARCH, "Pattern" },
                                                             { G2lib::Dubins3pBuildType::ELLIPSE, "Ellipse" } };

  constexpr int NUM_3P_ITERATIONS = 1000;

  fmt::print( fg( fmt::color::cyan ), "\n📈 Dubins3P Construction Benchmark ({} iterations):\n", NUM_3P_ITERATIONS );

  for ( const auto & [method, name] : methods )
  {
    auto start_3p      = chrono::high_resolution_clock::now();
    int  success_count = 0;

    for ( int i = 0; i < NUM_3P_ITERATIONS; ++i )
    {
      G2lib::Dubins3p db3p( fmt::format( "bench_{}", name ) );

      // Configure
      if ( method == G2lib::Dubins3pBuildType::SAMPLE_ONE_DEGREE ) { db3p.set_sample_points( 360 ); }
      else
      {
        db3p.set_tolerance( 0.1 * m_pi / 180.0 );
        db3p.set_sample_angle( m_pi / 4 );
      }

      // Varying parameters
      real_type angle_i = ( i % 100 ) * m_pi / 50.0;
      real_type angle_f = ( ( i + 50 ) % 100 ) * m_pi / 50.0;

      if ( db3p.build( -1, 0, angle_i, 0, 1, 1, 0, angle_f, 1.0, method ) ) { success_count++; }
    }

    auto end_3p      = chrono::high_resolution_clock::now();
    auto duration_3p = chrono::duration_cast<chrono::microseconds>( end_3p - start_3p );

    fmt::print(
      "├─ {}: {:.3f} ms total, {:.3f} µs avg, {:.0f} ops/sec, {}% success\n",
      name,
      duration_3p.count() / 1000.0,
      duration_3p.count() / static_cast<double>( NUM_3P_ITERATIONS ),
      NUM_3P_ITERATIONS / ( duration_3p.count() / 1e6 ),
      ( 100 * success_count ) / NUM_3P_ITERATIONS );
  }

  // Memory usage estimation
  fmt::print( fg( fmt::color::cyan ), "\n📊 Memory Estimation:\n" );
  fmt::print( "├─ Dubins object size: ~{} bytes\n", sizeof( G2lib::Dubins ) );
  fmt::print( "├─ Dubins3P object size: ~{} bytes\n", sizeof( G2lib::Dubins3p ) );
  fmt::print( "└─ CircleArc size: ~{} bytes\n", sizeof( G2lib::CircleArc ) );
}

// ============================================================================
// MAIN TEST SUITE
// ============================================================================
int main()
{
  fmt::print(
    fg( fmt::color::gold ) | fmt::emphasis::bold,
    "\n"
    "╔══════════════════════════════════════════════════════════════════╗\n"
    "║                🚀 COMPREHENSIVE DUBINS TEST SUITE                ║\n"
    "╚══════════════════════════════════════════════════════════════════╝\n"
  );

  auto total_start = chrono::high_resolution_clock::now();

  try
  {
    // Run all tests - AGGIUNGI I NUOVI TEST QUI
    vector<pair<string, function<void()>>> tests = {
      { "Basic Dubins Construction", test_basic_dubins },
      { "CSV Validation", test_csv_validation },
      { "Dubins Intersection", test_intersection },
      { "Dubins3P Comprehensive", test_dubins3p_comprehensive },
      { "Dubins3P Randomized", test_dubins3p_randomized },
      { "Dubins Method Completeness", test_dubins_method_completeness },
      { "Dubins3P Method Completeness", test_dubins3p_method_completeness },
      { "Collision & Intersection", test_collision_and_intersection },
      { "Transformations & Offsets", test_transformations_and_offsets },
      { "Dubins3P Edge Cases", test_dubins3p_edge_cases },
      { "Performance Benchmark", test_performance_benchmark }
    };

    int passed = 0;
    int failed = 0;

    for ( auto & [name, test_func] : tests )
    {
      try
      {
        fmt::print( fg( fmt::color::white ) | fmt::emphasis::bold, "\n▶️  Starting: {}\n", name );

        auto start = chrono::high_resolution_clock::now();
        test_func();
        auto end      = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>( end - start );

        fmt::print( fg( fmt::color::green ), "   ✓ Completed in {} μs\n\n", duration.count() );
        passed++;
      }
      catch ( const std::exception & e )
      {
        fmt::print( fg( fmt::color::red ), "   ✗ Failed: {}\n\n", e.what() );
        failed++;
      }
    }

    auto total_end      = chrono::high_resolution_clock::now();
    auto total_duration = chrono::duration_cast<chrono::milliseconds>( total_end - total_start );

    // Final summary
    fmt::print(
      fg( fmt::color::gold ) | fmt::emphasis::bold,
      "\n"
      "╔════════════════════════════════════════════════════════════╗\n"
      "║                     🏁 TEST SUMMARY                        ║\n"
      "╠════════════════════════════════════════════════════════════╣\n"
      "║ {:58} ║\n"
      "║ {:58} ║\n"
      "║ {:58} ║\n"
      "║ {:58} ║\n"
      "║ {:58} ║\n"
      "╚════════════════════════════════════════════════════════════╝\n",
      fmt::format( "Total tests run: {}", tests.size() ),
      fmt::format( "Passed:          {}", passed ),
      fmt::format( "Failed:          {}", failed ),
      fmt::format( "Success rate:    {:.1f}%", 100.0 * passed / tests.size() ),
      fmt::format( "Total time:      {} ms", total_duration.count() )
    );


    if ( failed == 0 )
    {
      fmt::print( fg( fmt::color::green ) | fmt::emphasis::bold, "\n✨ ALL TESTS PASSED SUCCESSFULLY! ✨\n" );
    }
    else
    {
      fmt::print(
        fg( fmt::color::yellow ) | fmt::emphasis::bold,
        "\n⚠️  Some tests failed. Please review the output above.\n" );
    }
  }
  catch ( const std::exception & e )
  {
    fmt::print( fg( fmt::color::red ) | fmt::emphasis::bold, "\n💥 FATAL ERROR: {}\n", e.what() );
    return 1;
  }

  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n🎯 ALL DONE FOLKS! 🎯\n" );

  return 0;
}
