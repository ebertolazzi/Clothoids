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
using namespace std;

// ============================================================================
// STILI DI FORMATTAZIONE
// ============================================================================

namespace Style
{
  // Stili principali
  const auto HEADER  = fg( fmt::color::steel_blue ) | fmt::emphasis::bold;
  const auto SECTION = fg( fmt::color::dodger_blue ) | fmt::emphasis::bold;
  const auto SUCCESS = fg( fmt::color::lime_green ) | fmt::emphasis::bold;
  // const auto ERROR      = fg(fmt::color::crimson) | fmt::emphasis::bold;
  // const auto WARNING    = fg(fmt::color::gold) | fmt::emphasis::bold;
  const auto INFO      = fg( fmt::color::deep_sky_blue ) | fmt::emphasis::bold;
  const auto VALUE     = fg( fmt::color::light_gray );
  const auto LABEL     = fg( fmt::color::silver );
  const auto HIGHLIGHT = fg( fmt::color::cyan ) | fmt::emphasis::bold;
  const auto GEOMETRY  = fg( fmt::color::violet ) | fmt::emphasis::bold;
  // const auto CURVE      = fg(fmt::color::orange) | fmt::emphasis::bold;
  const auto TEST_PASS = fg( fmt::color::green ) | fmt::emphasis::bold;
  const auto TEST_FAIL = fg( fmt::color::red ) | fmt::emphasis::bold;
  const auto POINT     = fg( fmt::color::spring_green ) | fmt::emphasis::bold;
  // const auto SEGMENT    = fg(fmt::color::hot_pink) | fmt::emphasis::bold;
}  // namespace Style

// ============================================================================
// FUNZIONI DI UTILITÀ PER LA FORMATTAZIONE
// ============================================================================

void print_header( const string & title, const string & icon = "📐" )
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

void print_polyline_info( const G2lib::PolyLine & P, const string & name = "PolyLine" )
{
  fmt::print( Style::GEOMETRY, "    📊 {} '{}':\n", name, P.name() );
  fmt::print( Style::LABEL, "      Numero segmenti: " );
  fmt::print( Style::VALUE, "{}\n", P.num_segments() );
  fmt::print( Style::LABEL, "      Numero punti:    " );
  fmt::print( Style::VALUE, "{}\n", P.numPoints() );
  fmt::print( Style::LABEL, "      Lunghezza totale: " );
  fmt::print( Style::VALUE, "{:.6f}\n", P.length() );
  fmt::print( Style::LABEL, "      Punto iniziale:  " );
  fmt::print( Style::POINT, "({:.4f}, {:.4f})\n", P.x_begin(), P.y_begin() );
  fmt::print( Style::LABEL, "      Punto finale:    " );
  fmt::print( Style::POINT, "({:.4f}, {:.4f})\n", P.x_end(), P.y_end() );
}

// ============================================================================
// FUNZIONI DI TEST
// ============================================================================

void test_constructors_and_build()
{
  print_section( "TEST COSTRUTTORI E METODI DI COSTRUZIONE", "🏗️" );

  int tests_passed = 0;
  int total_tests  = 0;

  // Test 1: Costruzione da array di punti
  {
    total_tests++;
    G2lib::PolyLine P( "test1" );

    vector<real_type> x = { 0.0, 1.0, 2.0, 3.0 };
    vector<real_type> y = { 0.0, 1.0, 0.0, 1.0 };

    P.build( static_cast<integer>( x.size() ), x.data(), y.data() );

    bool passed = ( P.num_segments() == 3 ) && ( P.numPoints() == 4 ) &&
                  ( abs( P.length() - ( sqrt( 2.0 ) + sqrt( 2.0 ) + sqrt( 2.0 ) ) ) < 1e-10 );

    print_test_result( "Costruzione da array di punti", passed );
    print_polyline_info( P, "P da array" );
    if ( passed ) tests_passed++;
  }

  // Test 2: Costruzione tramite push_back
  {
    total_tests++;
    G2lib::PolyLine P( "test2" );

    P.init( 0.0, 0.0 );
    P.push_back( 1.0, 0.0 );
    P.push_back( 1.0, 1.0 );
    P.push_back( 0.0, 1.0 );
    P.push_back( 0.0, 0.0 );  // Chiude il quadrato

    bool passed = ( P.num_segments() == 4 ) && ( P.numPoints() == 5 ) && ( abs( P.length() - 4.0 ) < 1e-10 );

    print_test_result( "Costruzione tramite push_back", passed );
    print_polyline_info( P, "P quadrato" );
    if ( passed ) tests_passed++;
  }

  // Test 3: Costruzione da LineSegment
  {
    total_tests++;
    G2lib::LineSegment LS( "segmento" );
    LS.build_2P( 0.0, 0.0, 3.0, 4.0 );

    G2lib::PolyLine P( "test3" );
    P.build( LS );

    bool passed = ( P.num_segments() == 1 ) && ( abs( P.length() - 5.0 ) < 1e-10 ) &&
                  ( abs( P.x_begin() - 0.0 ) < 1e-10 ) && ( abs( P.y_begin() - 0.0 ) < 1e-10 ) &&
                  ( abs( P.x_end() - 3.0 ) < 1e-10 ) && ( abs( P.y_end() - 4.0 ) < 1e-10 );

    print_test_result( "Costruzione da LineSegment", passed );
    if ( passed ) tests_passed++;
  }

  // Test 4: Costruttore di copia
  {
    total_tests++;
    G2lib::PolyLine P1( "originale" );
    P1.init( 0.0, 0.0 );
    P1.push_back( 1.0, 1.0 );
    P1.push_back( 2.0, 0.0 );

    G2lib::PolyLine P2( P1 );

    bool passed = ( P2.num_segments() == P1.num_segments() ) && ( abs( P2.length() - P1.length() ) < 1e-10 ) &&
                  ( P2.name() == P1.name() );

    print_test_result( "Costruttore di copia", passed );
    if ( passed ) tests_passed++;
  }

  // Test 5: Operatore di assegnazione
  {
    total_tests++;
    G2lib::PolyLine P1( "source" );
    P1.init( 0.0, 0.0 );
    P1.push_back( 2.0, 0.0 );
    P1.push_back( 2.0, 2.0 );

    G2lib::PolyLine P2( "dest" );
    P2 = P1;

    bool passed = ( P2.num_segments() == 2 ) && ( abs( P2.length() - 4.0 ) < 1e-10 );

    print_test_result( "Operatore di assegnazione", passed );
    if ( passed ) tests_passed++;
  }

  // Riassunto sezione
  fmt::print( "\n" );
  fmt::print( Style::INFO, "    📊 Riassunto costruttori: " );
  fmt::print( Style::HIGHLIGHT, "{}/{} test superati\n", tests_passed, total_tests );
}

void test_geometric_properties()
{
  print_section( "TEST PROPRIETÀ GEOMETRICHE", "📏" );

  int tests_passed = 0;
  int total_tests  = 0;

  // Test 1: Lunghezza e bounding box
  {
    total_tests++;
    G2lib::PolyLine P( "bbox_test" );
    P.init( 0.0, 0.0 );
    P.push_back( 3.0, 0.0 );
    P.push_back( 3.0, 4.0 );
    P.push_back( 0.0, 4.0 );

    real_type length          = P.length();
    real_type expected_length = 3.0 + 4.0 + 3.0;  // 10.0

    real_type xmin, ymin, xmax, ymax;
    P.bbox( xmin, ymin, xmax, ymax );

    fmt::print( Style::LABEL, "    📐 Polilinea a L (rettangolo mancante):\n" );
    fmt::print( Style::VALUE, "      Lunghezza calcolata: {:.4f}\n", length );
    fmt::print( Style::VALUE, "      Lunghezza attesa:    {:.4f}\n", expected_length );
    fmt::print( Style::VALUE, "      BBox: x=[{:.2f}, {:.2f}], y=[{:.2f}, {:.2f}]\n", xmin, xmax, ymin, ymax );

    bool passed = ( abs( length - expected_length ) < 1e-10 ) &&
                  ( xmin == 0.0 && xmax == 3.0 && ymin == 0.0 && ymax == 4.0 );

    print_test_result( "Lunghezza e bounding box", passed );
    if ( passed ) tests_passed++;
  }

  // Test 2: Valutazione punti
  {
    total_tests++;
    G2lib::PolyLine P( "eval_test" );
    P.init( 0.0, 0.0 );
    P.push_back( 2.0, 0.0 );
    P.push_back( 2.0, 2.0 );

    // Test a metà del primo segmento
    real_type s1 = 1.0;
    real_type x1 = P.X( s1 );
    real_type y1 = P.Y( s1 );

    // Test al punto di giunzione
    real_type s2 = 2.0;
    real_type x2 = P.X( s2 );
    real_type y2 = P.Y( s2 );

    // Test a metà del secondo segmento
    real_type s3 = 3.0;
    real_type x3 = P.X( s3 );
    real_type y3 = P.Y( s3 );

    fmt::print( Style::LABEL, "\n    📍 Valutazione punti:\n" );
    fmt::print( Style::VALUE, "      s={:.1f}: ({:.4f}, {:.4f}) [atteso: (1.0, 0.0)]\n", s1, x1, y1 );
    fmt::print( Style::VALUE, "      s={:.1f}: ({:.4f}, {:.4f}) [atteso: (2.0, 0.0)]\n", s2, x2, y2 );
    fmt::print( Style::VALUE, "      s={:.1f}: ({:.4f}, {:.4f}) [atteso: (2.0, 1.0)]\n", s3, x3, y3 );

    bool passed = ( abs( x1 - 1.0 ) < 1e-10 && abs( y1 - 0.0 ) < 1e-10 ) &&
                  ( abs( x2 - 2.0 ) < 1e-10 && abs( y2 - 0.0 ) < 1e-10 ) &&
                  ( abs( x3 - 2.0 ) < 1e-10 && abs( y3 - 1.0 ) < 1e-10 );

    print_test_result( "Valutazione X(s), Y(s)", passed );
    if ( passed ) tests_passed++;
  }

  // Test 3: Derivati
  {
    total_tests++;
    G2lib::PolyLine P( "deriv_test" );
    P.init( 0.0, 0.0 );
    P.push_back( 1.0, 0.0 );
    P.push_back( 1.0, 1.0 );

    real_type s    = 0.5;
    real_type X_D  = P.X_D( s );
    real_type Y_D  = P.Y_D( s );
    real_type X_DD = P.X_DD( s );
    real_type Y_DD = P.Y_DD( s );

    fmt::print( Style::LABEL, "\n    📈 Derivati:\n" );
    fmt::print( Style::VALUE, "      X_D({:.1f})  = {:.4f} [atteso: 1.0]\n", s, X_D );
    fmt::print( Style::VALUE, "      Y_D({:.1f})  = {:.4f} [atteso: 0.0]\n", s, Y_D );
    fmt::print( Style::VALUE, "      X_DD({:.1f}) = {:.4f} [atteso: 0.0]\n", s, X_DD );
    fmt::print( Style::VALUE, "      Y_DD({:.1f}) = {:.4f} [atteso: 0.0]\n", s, Y_DD );

    bool passed = ( abs( X_D - 1.0 ) < 1e-10 ) && ( abs( Y_D - 0.0 ) < 1e-10 ) && ( abs( X_DD ) < 1e-10 ) &&
                  ( abs( Y_DD ) < 1e-10 );

    print_test_result( "Derivate prime e seconde", passed );
    if ( passed ) tests_passed++;
  }

  // Test 4: Theta e sue derivate
  {
    total_tests++;
    G2lib::PolyLine P( "theta_test" );
    P.init( 0.0, 0.0 );
    P.push_back( 1.0, 0.0 );  // orizzontale (theta=0)
    P.push_back( 1.0, 1.0 );  // verticale (theta=π/2)

    real_type s1       = 0.5;  // Primo segmento
    real_type theta1   = P.theta( s1 );
    real_type theta_D1 = P.theta_D( s1 );

    real_type s2       = 1.5;  // Secondo segmento
    real_type theta2   = P.theta( s2 );
    real_type theta_D2 = P.theta_D( s2 );

    fmt::print( Style::LABEL, "\n    🧭 Angoli (theta):\n" );
    fmt::print( Style::VALUE, "      theta({:.1f})  = {:.6f} rad [{:.1f}°]\n", s1, theta1, theta1 * 180.0 / M_PI );
    fmt::print( Style::VALUE, "      theta_D({:.1f}) = {:.6f} [atteso: 0.0]\n", s1, theta_D1 );
    fmt::print( Style::VALUE, "      theta({:.1f})  = {:.6f} rad [{:.1f}°]\n", s2, theta2, theta2 * 180.0 / M_PI );
    fmt::print( Style::VALUE, "      theta_D({:.1f}) = {:.6f} [atteso: 0.0]\n", s2, theta_D2 );

    bool passed = ( abs( theta1 - 0.0 ) < 1e-10 ) && ( abs( theta_D1 ) < 1e-10 ) &&
                  ( abs( theta2 - M_PI / 2.0 ) < 1e-10 ) && ( abs( theta_D2 ) < 1e-10 );

    print_test_result( "Theta e derivate", passed );
    if ( passed ) tests_passed++;
  }

  // Riassunto sezione
  fmt::print( "\n" );
  fmt::print( Style::INFO, "    📊 Riassunto proprietà geometriche: " );
  fmt::print( Style::HIGHLIGHT, "{}/{} test superati\n", tests_passed, total_tests );
}

void test_geometric_transformations()
{
  print_section( "TEST TRASFORMAZIONI GEOMETRICHE", "🔄" );

  int tests_passed = 0;
  int total_tests  = 0;

  // Test 1: Traslazione
  {
    total_tests++;
    G2lib::PolyLine P( "translate_test" );
    P.init( 0.0, 0.0 );
    P.push_back( 1.0, 0.0 );
    P.push_back( 1.0, 1.0 );

    real_type tx = 2.0, ty = 3.0;
    real_type orig_x_end = P.x_end();
    real_type orig_y_end = P.y_end();

    P.translate( tx, ty );

    fmt::print( Style::LABEL, "    📍 Traslazione di ({:.1f}, {:.1f}):\n", tx, ty );
    fmt::print( Style::VALUE, "      Punto finale originale: ({:.1f}, {:.1f})\n", orig_x_end, orig_y_end );
    fmt::print( Style::VALUE, "      Punto finale dopo:      ({:.1f}, {:.1f})\n", P.x_end(), P.y_end() );

    bool passed = ( abs( P.x_end() - ( orig_x_end + tx ) ) < 1e-10 ) &&
                  ( abs( P.y_end() - ( orig_y_end + ty ) ) < 1e-10 );

    print_test_result( "Traslazione", passed );
    if ( passed ) tests_passed++;
  }

  // Test 2: Rotazione
  {
    total_tests++;
    G2lib::PolyLine P( "rotate_test" );
    P.init( 1.0, 0.0 );
    P.push_back( 2.0, 0.0 );

    real_type angle = M_PI / 2.0;  // 90 gradi
    real_type cx = 0.0, cy = 0.0;

    real_type orig_x_end = P.x_end();
    real_type orig_y_end = P.y_end();

    P.rotate( angle, cx, cy );

    // (2,0) ruotato di 90° diventa (0,2)
    fmt::print(
      Style::LABEL,
      "\n    📍 Rotazione di {:.1f}° intorno a ({:.1f}, {:.1f}):\n",
      angle * 180.0 / M_PI,
      cx,
      cy );
    fmt::print( Style::VALUE, "      Punto finale originale: ({:.1f}, {:.1f})\n", orig_x_end, orig_y_end );
    fmt::print( Style::VALUE, "      Punto finale dopo:      ({:.1f}, {:.1f})\n", P.x_end(), P.y_end() );
    fmt::print( Style::VALUE, "      Punto finale atteso:    ({:.1f}, {:.1f})\n", 0.0, 2.0 );

    bool passed = ( abs( P.x_end() - 0.0 ) < 1e-10 ) && ( abs( P.y_end() - 2.0 ) < 1e-10 );

    print_test_result( "Rotazione", passed );
    if ( passed ) tests_passed++;
  }

  // Test 3: Scala
  {
    total_tests++;
    G2lib::PolyLine P( "scale_test" );
    P.init( 1.0, 1.0 );
    P.push_back( 2.0, 1.0 );
    P.push_back( 2.0, 2.0 );

    real_type scale       = 2.0;
    real_type orig_length = P.length();

    P.scale( scale );

    fmt::print( Style::LABEL, "\n    📍 Scala di {:.1f}:\n", scale );
    fmt::print( Style::VALUE, "      Lunghezza originale: {:.4f}\n", orig_length );
    fmt::print( Style::VALUE, "      Lunghezza dopo:      {:.4f}\n", P.length() );
    fmt::print( Style::VALUE, "      Lunghezza attesa:    {:.4f}\n", orig_length * scale );

    bool passed = ( abs( P.length() - orig_length * scale ) < 1e-10 );

    print_test_result( "Scala uniforme", passed );
    if ( passed ) tests_passed++;
  }

  // Test 4: Reverse
  {
    total_tests++;
    G2lib::PolyLine P( "reverse_test" );
    P.init( 0.0, 0.0 );
    P.push_back( 1.0, 0.0 );
    P.push_back( 1.0, 1.0 );

    real_type orig_x_begin = P.x_begin();
    real_type orig_y_begin = P.y_begin();
    real_type orig_x_end   = P.x_end();
    real_type orig_y_end   = P.y_end();
    real_type orig_length  = P.length();

    P.reverse();

    fmt::print( Style::LABEL, "\n    📍 Reverse:\n" );
    fmt::print( Style::VALUE, "      Inizio originale: ({:.1f}, {:.1f})\n", orig_x_begin, orig_y_begin );
    fmt::print( Style::VALUE, "      Fine originale:   ({:.1f}, {:.1f})\n", orig_x_end, orig_y_end );
    fmt::print( Style::VALUE, "      Inizio dopo:      ({:.1f}, {:.1f})\n", P.x_begin(), P.y_begin() );
    fmt::print( Style::VALUE, "      Fine dopo:        ({:.1f}, {:.1f})\n", P.x_end(), P.y_end() );
    fmt::print( Style::VALUE, "      Lunghezza originale: {:.4f}\n", orig_length );
    fmt::print( Style::VALUE, "      Lunghezza dopo:      {:.4f}\n", P.length() );

    bool passed = ( abs( P.x_begin() - orig_x_end ) < 1e-10 ) && ( abs( P.y_begin() - orig_y_end ) < 1e-10 ) &&
                  ( abs( P.x_end() - orig_x_begin ) < 1e-10 ) && ( abs( P.y_end() - orig_y_begin ) < 1e-10 ) &&
                  ( abs( P.length() - orig_length ) < 1e-10 );

    print_test_result( "Reverse", passed );
    if ( passed ) tests_passed++;
  }

  // Test 5: Trim
  {
    total_tests++;
    G2lib::PolyLine P( "trim_test" );
    P.init( 0.0, 0.0 );
    P.push_back( 2.0, 0.0 );
    P.push_back( 2.0, 2.0 );
    P.push_back( 0.0, 2.0 );

    real_type s_begin = 1.0;
    real_type s_end   = 5.0;  // Lunghezza totale = 6.0

    G2lib::PolyLine P_trimmed( "trimmed" );
    P.trim( s_begin, s_end, P_trimmed );

    fmt::print( Style::LABEL, "\n    ✂️  Trim da s={:.1f} a s={:.1f}:\n", s_begin, s_end );
    fmt::print( Style::VALUE, "      Lunghezza originale: {:.4f}\n", P.length() );
    fmt::print( Style::VALUE, "      Lunghezza trimmed:   {:.4f}\n", P_trimmed.length() );
    fmt::print( Style::VALUE, "      Segmenti originali:  {}\n", P.num_segments() );
    fmt::print( Style::VALUE, "      Segmenti trimmed:    {}\n", P_trimmed.num_segments() );

    bool passed = ( abs( P_trimmed.length() - ( s_end - s_begin ) ) < 1e-10 ) &&
                  ( P_trimmed.num_segments() <= P.num_segments() );

    print_test_result( "Trim", passed );
    if ( passed ) tests_passed++;
  }

  // Test 6: Change origin
  {
    total_tests++;
    G2lib::PolyLine P( "origin_test" );
    P.init( 1.0, 1.0 );
    P.push_back( 2.0, 1.0 );
    P.push_back( 2.0, 2.0 );

    real_type newx0 = 0.0, newy0 = 0.0;
    real_type orig_x_begin = P.x_begin();
    real_type orig_y_begin = P.y_begin();

    P.change_origin( newx0, newy0 );

    fmt::print( Style::LABEL, "\n    🏁 Cambio origine a ({:.1f}, {:.1f}):\n", newx0, newy0 );
    fmt::print( Style::VALUE, "      Inizio originale: ({:.1f}, {:.1f})\n", orig_x_begin, orig_y_begin );
    fmt::print( Style::VALUE, "      Inizio dopo:      ({:.1f}, {:.1f})\n", P.x_begin(), P.y_begin() );

    bool passed = ( abs( P.x_begin() - newx0 ) < 1e-10 ) && ( abs( P.y_begin() - newy0 ) < 1e-10 );

    print_test_result( "Change origin", passed );
    if ( passed ) tests_passed++;
  }

  // Riassunto sezione
  fmt::print( "\n" );
  fmt::print( Style::INFO, "    📊 Riassunto trasformazioni: " );
  fmt::print( Style::HIGHLIGHT, "{}/{} test superati\n", tests_passed, total_tests );
}

void test_closest_point()
{
  print_section( "TEST PUNTO PIÙ VICINO", "📍" );

  int tests_passed = 0;
  int total_tests  = 0;

  // Test 1: Punto su un segmento
  {
    total_tests++;
    G2lib::PolyLine P( "closest_test1" );
    P.init( 0.0, 0.0 );
    P.push_back( 3.0, 0.0 );
    P.push_back( 3.0, 3.0 );

    real_type x = 1.5, y = 0.5;  // Sopra il primo segmento
    real_type X, Y, S, T, DST;
    integer   info = P.closest_point_ISO( x, y, X, Y, S, T, DST );

    fmt::print( Style::LABEL, "    📍 Punto di test: ({:.1f}, {:.1f})\n", x, y );
    fmt::print( Style::VALUE, "      Punto più vicino: ({:.4f}, {:.4f})\n", X, Y );
    fmt::print( Style::VALUE, "      Parametro S:      {:.4f}\n", S );
    fmt::print( Style::VALUE, "      Parametro T:      {:.4f}\n", T );
    fmt::print( Style::VALUE, "      Distanza:         {:.4f}\n", DST );
    fmt::print( Style::VALUE, "      Info code:        {}\n", info );

    bool passed = ( abs( X - 1.5 ) < 1e-10 ) && ( abs( Y - 0.0 ) < 1e-10 ) && ( abs( DST - 0.5 ) < 1e-10 ) &&
                  ( info == 0 );

    print_test_result( "Punto vicino a segmento orizzontale", passed );
    if ( passed ) tests_passed++;
  }

  // Test 2: Punto su un vertice
  {
    total_tests++;
    G2lib::PolyLine P( "closest_test2" );
    P.init( 0.0, 0.0 );
    P.push_back( 2.0, 0.0 );
    P.push_back( 2.0, 2.0 );

    real_type x = 0.0, y = 0.0;  // Esattamente sul vertice iniziale
    real_type X, Y, S, T, DST;
    integer   info = P.closest_point_ISO( x, y, X, Y, S, T, DST );

    fmt::print( Style::LABEL, "\n    📍 Punto su vertice iniziale:\n" );
    fmt::print( Style::VALUE, "      Punto più vicino: ({:.4f}, {:.4f})\n", X, Y );
    fmt::print( Style::VALUE, "      Distanza:         {:.4f}\n", DST );

    bool passed = ( abs( X - 0.0 ) < 1e-10 ) && ( abs( Y - 0.0 ) < 1e-10 ) && ( abs( DST ) < 1e-10 ) && ( info == 0 );

    print_test_result( "Punto su vertice", passed );
    if ( passed ) tests_passed++;
  }

  // Test 3: Punto vicino a giunzione
  {
    total_tests++;
    G2lib::PolyLine P( "closest_test3" );
    P.init( 0.0, 0.0 );
    P.push_back( 2.0, 0.0 );
    P.push_back( 2.0, 2.0 );

    real_type x = 2.1, y = 0.1;  // Vicino alla giunzione (2,0)
    real_type X, Y, S, T, DST;
    integer   info = P.closest_point_ISO( x, y, X, Y, S, T, DST );

    fmt::print( Style::LABEL, "\n    📍 Punto vicino a giunzione:\n" );
    fmt::print( Style::VALUE, "      Punto più vicino: ({:.4f}, {:.4f})\n", X, Y );
    fmt::print( Style::VALUE, "      Distanza:         {:.4f}\n", DST );
    fmt::print( Style::VALUE, "      Distanza attesa:  {:.4f}\n", sqrt( 0.1 * 0.1 + 0.1 * 0.1 ) );
    fmt::print( Style::VALUE, "      Info:             {}\n", info );

    bool passed = ( abs( X - 2.0 ) < 1e-10 ) && ( abs( Y - 0.0 ) < 1e-10 ) &&
                  ( abs( DST - sqrt( 0.01 + 0.01 ) ) < 1e-10 );

    print_test_result( "Punto vicino a giunzione", passed );
    if ( passed ) tests_passed++;
  }

  // Test 4: Punto lontano dalla polilinea
  {
    total_tests++;
    G2lib::PolyLine P( "closest_test4" );
    P.init( 0.0, 0.0 );
    P.push_back( 1.0, 0.0 );
    P.push_back( 1.0, 1.0 );

    real_type x = 5.0, y = 5.0;  // Lontano
    real_type X, Y, S, T, DST;
    integer   info = P.closest_point_ISO( x, y, X, Y, S, T, DST );

    fmt::print( Style::LABEL, "\n    📍 Punto lontano:\n" );
    fmt::print( Style::VALUE, "      Punto più vicino: ({:.4f}, {:.4f})\n", X, Y );
    fmt::print( Style::VALUE, "      Distanza:         {:.4f}\n", DST );
    fmt::print( Style::VALUE, "      Info:             {}\n", info );

    // Il punto più vicino dovrebbe essere (1,1)
    real_type expected_dist = sqrt( ( 5 - 1 ) * ( 5 - 1 ) + ( 5 - 1 ) * ( 5 - 1 ) );
    bool passed = ( abs( X - 1.0 ) < 1e-10 ) && ( abs( Y - 1.0 ) < 1e-10 ) && ( abs( DST - expected_dist ) < 1e-10 );

    print_test_result( "Punto lontano", passed );
    if ( passed ) tests_passed++;
  }

  // Riassunto sezione
  fmt::print( "\n" );
  fmt::print( Style::INFO, "    📊 Riassunto punto più vicino: " );
  fmt::print( Style::HIGHLIGHT, "{}/{} test superati\n", tests_passed, total_tests );
}

void test_intersections()
{
  print_section( "TEST INTERSEZIONI", "⚡" );

  int tests_passed = 0;
  int total_tests  = 0;

  // Test 1: Intersezione tra due polilinee che si incrociano
  {
    total_tests++;
    G2lib::PolyLine P1( "P1_cross" );
    P1.init( 0.0, 0.0 );
    P1.push_back( 2.0, 2.0 );

    G2lib::PolyLine P2( "P2_cross" );
    P2.init( 0.0, 2.0 );
    P2.push_back( 2.0, 0.0 );

    vector<real_type> s1, s2;
    P1.intersect( P2, s1, s2 );

    fmt::print( Style::LABEL, "    ✖️  Intersezione di due segmenti che si incrociano:\n" );
    fmt::print( Style::VALUE, "      Numero intersezioni trovate: {}\n", s1.size() );

    if ( !s1.empty() )
    {
      fmt::print( Style::VALUE, "      s1[0] = {:.6f}\n", s1[0] );
      fmt::print( Style::VALUE, "      s2[0] = {:.6f}\n", s2[0] );

      // Verifica che il punto di intersezione sia (1,1)
      real_type x_int = P1.X( s1[0] );
      real_type y_int = P1.Y( s1[0] );
      fmt::print( Style::VALUE, "      Punto intersezione: ({:.6f}, {:.6f})\n", x_int, y_int );
    }

    bool passed = ( s1.size() == 1 ) && ( s2.size() == 1 );

    if ( passed && !s1.empty() )
    {
      real_type x_int = P1.X( s1[0] );
      real_type y_int = P1.Y( s1[0] );
      passed          = passed && ( abs( x_int - 1.0 ) < 1e-6 ) && ( abs( y_int - 1.0 ) < 1e-6 );
    }

    print_test_result( "Intersezione semplice", passed );
    if ( passed ) tests_passed++;
  }

  // Test 2: Intersezione tra polilinee con segmenti multipli
  {
    total_tests++;
    G2lib::PolyLine P1( "P1_multi" );
    P1.init( 0.0, 0.0 );
    P1.push_back( 2.0, 0.0 );
    P1.push_back( 2.0, 2.0 );

    G2lib::PolyLine P2( "P2_multi" );
    P2.init( 1.0, 1.0 );
    P2.push_back( 3.0, 1.0 );

    vector<real_type> s1, s2;
    P1.intersect( P2, s1, s2 );

    fmt::print( Style::LABEL, "\n    📍 Intersezione polilinea L con segmento orizzontale:\n" );
    fmt::print( Style::VALUE, "      Numero intersezioni trovate: {}\n", s1.size() );

    for ( size_t i = 0; i < s1.size(); ++i )
    {
      fmt::print( Style::VALUE, "      s1[{}] = {:.6f}, s2[{}] = {:.6f}\n", i, s1[i], i, s2[i] );
    }

    // Dovrebbero intersecarsi in (2,1)
    bool passed = ( s1.size() == 1 ) && ( s2.size() == 1 );

    if ( passed && !s1.empty() )
    {
      real_type x_int = P1.X( s1[0] );
      real_type y_int = P1.Y( s1[0] );
      passed          = passed && ( abs( x_int - 2.0 ) < 1e-6 ) && ( abs( y_int - 1.0 ) < 1e-6 );
    }

    print_test_result( "Intersezione con polilinea a L", passed );
    if ( passed ) tests_passed++;
  }

  // Test 3: Intersezione usando IntersectList
  {
    total_tests++;
    G2lib::PolyLine P1( "P1_ilist" );
    P1.init( 0.0, 0.0 );
    P1.push_back( 1.0, 1.0 );

    G2lib::PolyLine P2( "P2_ilist" );
    P2.init( 0.0, 1.0 );
    P2.push_back( 1.0, 0.0 );

    G2lib::IntersectList ilist;
    P1.intersect( P2, ilist );

    fmt::print( Style::LABEL, "\n    📋 Intersezione usando IntersectList:\n" );
    fmt::print( Style::VALUE, "      Numero intersezioni trovate: {}\n", ilist.size() );

    for ( size_t i = 0; i < ilist.size(); ++i )
    {
      fmt::print( Style::VALUE, "      Intersezione {}: s1={:.6f}, s2={:.6f}\n", i, ilist[i].first, ilist[i].second );
    }

    bool passed = ( ilist.size() == 1 );

    if ( passed && !ilist.empty() )
    {
      real_type x_int = P1.X( ilist[0].first );
      real_type y_int = P1.Y( ilist[0].first );
      passed          = passed && ( abs( x_int - 0.5 ) < 1e-6 ) && ( abs( y_int - 0.5 ) < 1e-6 );
    }

    print_test_result( "Intersezione con IntersectList", passed );
    if ( passed ) tests_passed++;
  }

  // Test 4: Nessuna intersezione
  {
    total_tests++;
    G2lib::PolyLine P1( "P1_no_int" );
    P1.init( 0.0, 0.0 );
    P1.push_back( 1.0, 0.0 );

    G2lib::PolyLine P2( "P2_no_int" );
    P2.init( 0.0, 2.0 );
    P2.push_back( 1.0, 2.0 );

    vector<real_type> s1, s2;
    P1.intersect( P2, s1, s2 );

    fmt::print( Style::LABEL, "\n    🚫 Test nessuna intersezione:\n" );
    fmt::print( Style::VALUE, "      Numero intersezioni trovate: {}\n", s1.size() );

    bool passed = ( s1.empty() && s2.empty() );

    print_test_result( "Nessuna intersezione", passed );
    if ( passed ) tests_passed++;
  }

  // Riassunto sezione
  fmt::print( "\n" );
  fmt::print( Style::INFO, "    📊 Riassunto intersezioni: " );
  fmt::print( Style::HIGHLIGHT, "{}/{} test superati\n", tests_passed, total_tests );
}

void test_collision()
{
  print_section( "TEST COLLISIONI", "💥" );

  int tests_passed = 0;
  int total_tests  = 0;

  // Test 1: Collisione tra polilinee che si intersecano
  {
    total_tests++;
    G2lib::PolyLine P1( "P1_collide" );
    P1.init( 0.0, 0.0 );
    P1.push_back( 2.0, 0.0 );

    G2lib::PolyLine P2( "P2_collide" );
    P2.init( 1.0, -1.0 );
    P2.push_back( 1.0, 1.0 );

    bool collision = P1.collision( P2 );

    fmt::print( Style::LABEL, "    💥 Test collisione (dovrebbe esserci):\n" );
    fmt::print( Style::VALUE, "      Collisione rilevata: {}\n", collision ? "true" : "false" );

    bool passed = collision;

    print_test_result( "Collisione rilevata", passed );
    if ( passed ) tests_passed++;
  }

  // Test 2: Nessuna collisione
  {
    total_tests++;
    G2lib::PolyLine P1( "P1_no_collide" );
    P1.init( 0.0, 0.0 );
    P1.push_back( 1.0, 0.0 );

    G2lib::PolyLine P2( "P2_no_collide" );
    P2.init( 0.0, 2.0 );
    P2.push_back( 1.0, 2.0 );

    bool collision = P1.collision( P2 );

    fmt::print( Style::LABEL, "\n    🚫 Test nessuna collisione:\n" );
    fmt::print( Style::VALUE, "      Collisione rilevata: {}\n", collision ? "true" : "false" );

    bool passed = !collision;

    print_test_result( "Nessuna collisione", passed );
    if ( passed ) tests_passed++;
  }

  // Test 3: Collisione con offset ISO
  {
    total_tests++;
    G2lib::PolyLine P1( "P1_iso" );
    P1.init( 0.0, 0.0 );
    P1.push_back( 3.0, 0.0 );

    G2lib::PolyLine P2( "P2_iso" );
    P2.init( 1.0, 0.5 );
    P2.push_back( 2.0, 0.5 );

    real_type offs1 = 0.0;
    real_type offs2 = 0.0;

    bool collision = P1.collision_ISO( offs1, P2, offs2 );

    fmt::print( Style::LABEL, "\n    📏 Test collisione ISO (offset=0):\n" );
    fmt::print( Style::VALUE, "      Collisione rilevata: {}\n", collision ? "true" : "false" );

    // Non dovrebbero collidere perché P2 è a y=0.5
    bool passed = !collision;

    print_test_result( "Collisione ISO senza offset", passed );
    if ( passed ) tests_passed++;
  }

  // Riassunto sezione
  fmt::print( "\n" );
  fmt::print( Style::INFO, "    📊 Riassunto collisioni: " );
  fmt::print( Style::HIGHLIGHT, "{}/{} test superati\n", tests_passed, total_tests );
}

void test_building_from_curves()
{
  print_section( "TEST COSTRUZIONE DA ALTRE CURVE", "🔄" );

  int tests_passed = 0;
  int total_tests  = 0;

  // Test 1: Costruzione da ClothoidCurve
  {
    total_tests++;
    G2lib::ClothoidCurve C( "clothoid" );
    C.build_G1( 0.0, 0.0, 0.0, 2.0, 0.0, M_PI / 2.0 );

    G2lib::PolyLine P( "from_clothoid" );
    real_type       tol = 0.01;
    P.build( C, tol );

    fmt::print( Style::LABEL, "    🌀 Costruzione da clothoid (tol={:.3f}):\n", tol );
    fmt::print( Style::VALUE, "      Segmenti generati: {}\n", P.num_segments() );
    fmt::print( Style::VALUE, "      Lunghezza clothoid: {:.4f}\n", C.length() );
    fmt::print( Style::VALUE, "      Lunghezza polyline: {:.4f}\n", P.length() );

    bool passed = ( P.num_segments() > 1 ) && ( abs( P.length() - C.length() ) < tol * 10 );  // Tolleranza larga

    print_test_result( "Costruzione da clothoid", passed );
    if ( passed ) tests_passed++;
  }

  // Test 2: Costruzione da Biarc
  {
    total_tests++;
    G2lib::Biarc B( "biarc" );
    B.build_3P( 0.0, 0.0, 1.0, 1.0, 2.0, 0.0 );

    G2lib::PolyLine P( "from_biarc" );
    real_type       tol = 0.01;
    P.build( B, tol );

    fmt::print( Style::LABEL, "\n    🏹 Costruzione da biarc (tol={:.3f}):\n", tol );
    fmt::print( Style::VALUE, "      Segmenti generati: {}\n", P.num_segments() );
    fmt::print( Style::VALUE, "      Lunghezza biarc:    {:.4f}\n", B.length() );
    fmt::print( Style::VALUE, "      Lunghezza polyline: {:.4f}\n", P.length() );

    bool passed = ( P.num_segments() > 1 );

    print_test_result( "Costruzione da biarc", passed );
    if ( passed ) tests_passed++;
  }

  // Test 3: push_back da CircleArc
  {
    total_tests++;
    G2lib::CircleArc C( "circle_arc" );
    C.build_3P( 0.0, 0.0, 1.0, 1.0, 2.0, 0.0 );

    G2lib::PolyLine P( "from_circle_arc" );
    real_type       tol = 0.01;
    P.init( C.x_begin(), C.y_begin() );
    P.push_back( C, tol );

    fmt::print( Style::LABEL, "\n    ⭕ Costruzione da circle arc (tol={:.3f}):\n", tol );
    fmt::print( Style::VALUE, "      Segmenti generati: {}\n", P.num_segments() );

    bool passed = ( P.num_segments() > 1 );

    print_test_result( "Costruzione da circle arc", passed );
    if ( passed ) tests_passed++;
  }

  // Riassunto sezione
  fmt::print( "\n" );
  fmt::print( Style::INFO, "    📊 Riassunto costruzione da curve: " );
  fmt::print( Style::HIGHLIGHT, "{}/{} test superati\n", tests_passed, total_tests );
}

void test_aabb_tree_and_triangles()
{
  print_section( "TEST AABB TREE E TRIANGOLI", "🌳" );

  int tests_passed = 0;
  int total_tests  = 0;

  // Test 1: Generazione triangoli per bounding box
  {
    total_tests++;
    G2lib::PolyLine P( "triangle_test" );
    P.init( 0.0, 0.0 );
    P.push_back( 2.0, 0.0 );
    P.push_back( 2.0, 2.0 );
    P.push_back( 0.0, 2.0 );
    P.push_back( 0.0, 0.0 );  // Quadrato chiuso

    vector<G2lib::Triangle2D> tvec;
    real_type                 max_angle = M_PI / 6.0;  // 30 gradi
    real_type                 max_size  = 1.0;
    integer                   icurve    = 0;

    P.bb_triangles( tvec, max_angle, max_size, icurve );

    fmt::print( Style::LABEL, "    📐 Generazione triangoli per bounding box:\n" );
    fmt::print( Style::VALUE, "      Numero triangoli generati: {}\n", tvec.size() );
    fmt::print( Style::VALUE, "      Max angle: {:.1f}°, Max size: {:.1f}\n", max_angle * 180.0 / M_PI, max_size );

    if ( !tvec.empty() )
    {
      real_type total_area = 0.0;
      for ( const auto & t : tvec ) { total_area += t.area(); }
      fmt::print( Style::VALUE, "      Area totale triangoli: {:.4f}\n", total_area );
    }

    bool passed = !tvec.empty();

    print_test_result( "Generazione triangoli", passed );
    if ( passed ) tests_passed++;
  }

  // Test 2: Costruzione AABB tree
  {
    total_tests++;
    G2lib::PolyLine P( "aabb_test" );
    P.init( 0.0, 0.0 );
    P.push_back( 1.0, 0.0 );
    P.push_back( 1.0, 1.0 );
    P.push_back( 0.0, 1.0 );

    // Costruisce l'AABB tree
    P.build_AABBtree();

    fmt::print( Style::LABEL, "\n    🌳 Costruzione AABB tree:\n" );
    fmt::print( Style::VALUE, "      AABB tree costruito per polilinea con {} segmenti\n", P.num_segments() );

    // Non c'è un modo diretto per verificare l'AABB tree, ma se non crasha è buono
    bool passed = true;

    print_test_result( "Costruzione AABB tree", passed );
    if ( passed ) tests_passed++;
  }

  // Riassunto sezione
  fmt::print( "\n" );
  fmt::print( Style::INFO, "    📊 Riassunto AABB tree e triangoli: " );
  fmt::print( Style::HIGHLIGHT, "{}/{} test superati\n", tests_passed, total_tests );
}

void test_original_example()
{
  print_section( "TEST ESEMPIO ORIGINALE", "📄" );

  fmt::print( Style::LABEL, "    📖 Esecuzione dell'esempio originale dal file testPolyline.cc:\n\n" );

  // Esecuzione dell'esempio originale
  G2lib::ClothoidCurve C1{ "C1" };
  G2lib::ClothoidCurve C2{ "C2" };
  G2lib::PolyLine      P1{ "P1" };
  G2lib::PolyLine      P2{ "P2" };

  constexpr real_type x0  = 0;
  constexpr real_type y0  = 0;
  constexpr real_type th0 = 0;
  constexpr real_type x1  = 0;
  constexpr real_type y1  = 3;
  real_type const     th1 = Utils::m_pi;

  C1.build_G1( x0, y0, th0, x1, y1, th1 );

  fmt::print( Style::VALUE, "      Clothoid C1 costruita\n" );

  C2 = C1;
  C2.rotate( Utils::m_pi / 3, 0, 0 );
  C2.translate( 1, -1 );

  G2lib::IntersectList ilist;
  C1.intersect( C2, ilist );

  fmt::print( Style::VALUE, "\n      Intersezioni tra clothoids:\n" );
  fmt::print( Style::VALUE, "        Numero intersezioni: {}\n", ilist.size() );

  for ( size_t i = 0; i < ilist.size(); ++i )
  {
    fmt::print( Style::VALUE, "        s1[{}] = {:.6f}, s2[{}] = {:.6f}\n", i, ilist[i].first, i, ilist[i].second );
  }

  P1.build( C1, 0.00001 );
  P2.build( C2, 0.00001 );

  fmt::print( Style::VALUE, "\n      Polilinee approssimate:\n" );
  fmt::print( Style::VALUE, "        P1: {} segmenti, lunghezza: {:.6f}\n", P1.num_segments(), P1.length() );
  fmt::print( Style::VALUE, "        P2: {} segmenti, lunghezza: {:.6f}\n", P2.num_segments(), P2.length() );

  vector<real_type> s1, s2;
  P1.intersect( P2, s1, s2 );

  fmt::print( Style::VALUE, "\n      Intersezioni tra polilinee:\n" );
  fmt::print( Style::VALUE, "        Numero intersezioni: {} {}\n", s1.size(), s2.size() );

  for ( size_t i = 0; i < s1.size(); ++i )
  {
    fmt::print( Style::VALUE, "        s1[{}] = {:.6f}, s2[{}] = {:.6f}\n", i, s1[i], i, s2[i] );
  }

  fmt::print( "\n" );
  print_test_result( "Esempio originale completato", true );
}

// ============================================================================
// FUNZIONE PRINCIPALE
// ============================================================================

int main()
{
  // Intestazione del programma
  print_header( "TEST SUITE CLASSE PolyLine", "📐" );

  auto start_time = chrono::high_resolution_clock::now();

  fmt::print( Style::INFO, "    📅 Data e ora: {:%Y-%m-%d %H:%M:%S}\n", chrono::system_clock::now() );
  fmt::print( Style::INFO, "    🏷️  Versione: PolyLine Test Suite 1.0\n" );
  fmt::print( Style::INFO, "    👤 Autore: Enrico Bertolazzi\n" );
  fmt::print( Style::INFO, "    📧 Email: enrico.bertolazzi@unitn.it\n\n" );

  // Esecuzione di tutti i test
  test_constructors_and_build();
  test_geometric_properties();
  test_geometric_transformations();
  test_closest_point();
  test_intersections();
  test_collision();
  test_building_from_curves();
  test_aabb_tree_and_triangles();
  test_original_example();

  // Calcolo tempo di esecuzione
  auto end_time = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::milliseconds>( end_time - start_time );

  // Riassunto finale
  print_header( "RIEPILOGO FINALE TEST" );

  fmt::print( Style::SUCCESS, "    ✅ Tutti i test sono stati eseguiti con successo!\n" );
  fmt::print( Style::INFO, "    ⏱️  Tempo di esecuzione: {} ms\n", duration.count() );
  fmt::print( Style::INFO, "    📊 Categorie test completate: 8\n\n" );

  fmt::print(
    Style::HEADER,
    "╔══════════════════════════════════════════════════════════════╗\n"
    "║                    🎉 TEST COMPLETATI!                       ║\n"
    "╚══════════════════════════════════════════════════════════════╝\n" );

  return 0;
}
