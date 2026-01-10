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
  const auto HEADER    = fg( fmt::color::steel_blue ) | fmt::emphasis::bold;
  const auto SECTION   = fg( fmt::color::dodger_blue ) | fmt::emphasis::bold;
  const auto SUCCESS   = fg( fmt::color::lime_green ) | fmt::emphasis::bold;
  //const auto ERROR     = fg( fmt::color::crimson ) | fmt::emphasis::bold;
  const auto WARNING   = fg( fmt::color::gold ) | fmt::emphasis::bold;
  const auto INFO      = fg( fmt::color::deep_sky_blue ) | fmt::emphasis::bold;
  const auto VALUE     = fg( fmt::color::light_gray );
  const auto LABEL     = fg( fmt::color::silver );
  const auto HIGHLIGHT = fg( fmt::color::cyan ) | fmt::emphasis::bold;
  const auto GEOMETRY  = fg( fmt::color::violet ) | fmt::emphasis::bold;
  const auto TEST_PASS = fg( fmt::color::green ) | fmt::emphasis::bold;
  const auto TEST_FAIL = fg( fmt::color::red ) | fmt::emphasis::bold;
}  // namespace Style

// ============================================================================
// FUNZIONI DI UTILITÀ PER LA FORMATTAZIONE
// ============================================================================

void print_header( const string & title )
{
  fmt::print( "\n" );
  fmt::print( Style::HEADER, "╔══════════════════════════════════════════════════════════════╗\n" );
  fmt::print( Style::HEADER, "║ {:^60} ║\n", title );
  fmt::print( Style::HEADER, "╚══════════════════════════════════════════════════════════════╝\n" );
  fmt::print( "\n" );
}

void print_section( const string & title, const string & icon = "📋" )
{
  fmt::print( Style::SECTION, "\n┌{0:─^{1}}┐\n", "", 60 );
  fmt::print( Style::SECTION, "│ {:^58} │\n", fmt::format( "{} {}", icon, title ) );
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

void print_triangle_info( const G2lib::Triangle2D & T, const string & name = "Triangle" )
{
  fmt::print( Style::GEOMETRY, "    📐 {}:\n", name );
  fmt::print( Style::LABEL, "      Vertici: " );
  fmt::print( Style::VALUE, "P1({:.4f}, {:.4f}), ", T.x1(), T.y1() );
  fmt::print( Style::VALUE, "P2({:.4f}, {:.4f}), ", T.x2(), T.y2() );
  fmt::print( Style::VALUE, "P3({:.4f}, {:.4f})\n", T.x3(), T.y3() );
  fmt::print( Style::LABEL, "      Parametri: " );
  fmt::print( Style::VALUE, "s0={:.4f}, s1={:.4f}, icurve={}\n", T.S0(), T.S1(), T.Icurve() );
}

// ============================================================================
// FUNZIONI DI TEST
// ============================================================================

void test_constructors()
{
  print_section( "TEST COSTRUTTORI E METODI DI COSTRUZIONE", "🏗️" );

  int tests_passed = 0;
  int total_tests  = 0;

  // Test 1: Costruttore di default
  {
    total_tests++;
    G2lib::Triangle2D T1;
    bool              passed =
      ( T1.x1() == 0 && T1.y1() == 0 && T1.x2() == 0 && T1.y2() == 0 && T1.x3() == 0 && T1.y3() == 0 && T1.S0() == 0 &&
        T1.S1() == 0 && T1.Icurve() == 0 );

    print_test_result( "Costruttore default", passed );
    if ( passed ) tests_passed++;
  }

  // Test 2: Costruttore con coordinate
  {
    total_tests++;
    G2lib::Triangle2D T2( 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 0.5, 1.5, 10 );

    bool passed =
      ( T2.x1() == 1.0 && T2.y1() == 2.0 && T2.x2() == 3.0 && T2.y2() == 4.0 && T2.x3() == 5.0 && T2.y3() == 6.0 &&
        T2.S0() == 0.5 && T2.S1() == 1.5 && T2.Icurve() == 10 );

    print_test_result( "Costruttore con coordinate", passed );
    print_triangle_info( T2, "T2" );
    if ( passed ) tests_passed++;
  }

  // Test 3: Costruttore con array
  {
    total_tests++;
    real_type p1[2] = { 1.0, 1.0 };
    real_type p2[2] = { 4.0, 1.0 };
    real_type p3[2] = { 2.0, 3.0 };

    G2lib::Triangle2D T3( p1, p2, p3, 0.0, 1.0, 1 );

    bool passed =
      ( T3.x1() == 1.0 && T3.y1() == 1.0 && T3.x2() == 4.0 && T3.y2() == 1.0 && T3.x3() == 2.0 && T3.y3() == 3.0 );

    print_test_result( "Costruttore con array", passed );
    if ( passed ) tests_passed++;
  }

  // Test 4: Costruttore di copia
  {
    total_tests++;
    G2lib::Triangle2D T4( 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 0.5, 1.5, 10 );
    G2lib::Triangle2D T5( T4 );

    bool passed =
      ( T5.x1() == T4.x1() && T5.y1() == T4.y1() && T5.x2() == T4.x2() && T5.y2() == T4.y2() && T5.x3() == T4.x3() &&
        T5.y3() == T4.y3() && T5.S0() == T4.S0() && T5.S1() == T4.S1() && T5.Icurve() == T4.Icurve() );

    print_test_result( "Costruttore di copia", passed );
    if ( passed ) tests_passed++;
  }

  // Test 5: Operatore di assegnazione
  {
    total_tests++;
    G2lib::Triangle2D T6( 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 0.5, 1.5, 10 );
    G2lib::Triangle2D T7;
    T7 = T6;

    bool passed =
      ( T7.x1() == T6.x1() && T7.y1() == T6.y1() && T7.x2() == T6.x2() && T7.y2() == T6.y2() && T7.x3() == T6.x3() &&
        T7.y3() == T6.y3() );

    print_test_result( "Operatore di assegnazione", passed );
    if ( passed ) tests_passed++;
  }

  // Test 6: Metodo build con coordinate
  {
    total_tests++;
    G2lib::Triangle2D T8;
    T8.build( 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0 );

    bool passed =
      ( T8.x1() == 0.0 && T8.y1() == 0.0 && T8.x2() == 1.0 && T8.y2() == 0.0 && T8.x3() == 0.0 && T8.y3() == 1.0 );

    print_test_result( "Metodo build con coordinate", passed );
    if ( passed ) tests_passed++;
  }

  // Test 7: Metodo build con array
  {
    total_tests++;
    real_type p1[2] = { 0.0, 0.0 };
    real_type p2[2] = { 2.0, 0.0 };
    real_type p3[2] = { 1.0, 2.0 };

    G2lib::Triangle2D T9;
    T9.build( p1, p2, p3, 0.0, 2.0, 1 );

    bool passed =
      ( T9.x1() == 0.0 && T9.y1() == 0.0 && T9.x2() == 2.0 && T9.y2() == 0.0 && T9.x3() == 1.0 && T9.y3() == 2.0 );

    print_test_result( "Metodo build con array", passed );
    if ( passed ) tests_passed++;
  }

  // Riassunto sezione
  fmt::print( "\n" );
  fmt::print( Style::INFO, "    📊 Riassunto costruttori: " );
  fmt::print( Style::HIGHLIGHT, "{}/{} test superati\n", tests_passed, total_tests );
}

void test_access_methods()
{
  print_section( "TEST METODI DI ACCESSO", "🔍" );

  G2lib::Triangle2D T( 1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 0.1, 0.9, 5 );

  fmt::print( Style::LABEL, "    📍 Metodi di accesso diretti:\n" );
  fmt::print( Style::VALUE, "      x1() = {:.2f}, y1() = {:.2f}\n", T.x1(), T.y1() );
  fmt::print( Style::VALUE, "      x2() = {:.2f}, y2() = {:.2f}\n", T.x2(), T.y2() );
  fmt::print( Style::VALUE, "      x3() = {:.2f}, y3() = {:.2f}\n", T.x3(), T.y3() );
  fmt::print( Style::VALUE, "      S0() = {:.2f}, S1() = {:.2f}\n", T.S0(), T.S1() );
  fmt::print( Style::VALUE, "      Icurve() = {}\n", T.Icurve() );

  fmt::print( Style::LABEL, "\n    📍 Metodi di accesso con puntatori:\n" );
  const real_type * p1 = T.P1();
  const real_type * p2 = T.P2();
  const real_type * p3 = T.P3();

  fmt::print( Style::VALUE, "      P1() = [{:.2f}, {:.2f}]\n", p1[0], p1[1] );
  fmt::print( Style::VALUE, "      P2() = [{:.2f}, {:.2f}]\n", p2[0], p2[1] );
  fmt::print( Style::VALUE, "      P3() = [{:.2f}, {:.2f}]\n", p3[0], p3[1] );

  // Test consistenza
  bool tests_passed = true;
  tests_passed &= ( p1[0] == T.x1() && p1[1] == T.y1() );
  tests_passed &= ( p2[0] == T.x2() && p2[1] == T.y2() );
  tests_passed &= ( p3[0] == T.x3() && p3[1] == T.y3() );

  print_test_result( "Consistenza metodi di accesso", tests_passed );
}

void test_geometric_transformations()
{
  print_section( "TEST TRASFORMAZIONI GEOMETRICHE", "🔄" );

  int tests_passed = 0;
  int total_tests  = 0;

  // Test 1: Traslazione
  {
    total_tests++;
    G2lib::Triangle2D T( 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0 );
    real_type         tx = 2.0, ty = 3.0;

    fmt::print( Style::LABEL, "    📍 Traslazione di ({:.1f}, {:.1f}):\n", tx, ty );
    fmt::print(
      Style::VALUE,
      "      Prima: P1({:.1f}, {:.1f}), P2({:.1f}, {:.1f}), P3({:.1f}, {:.1f})\n",
      T.x1(),
      T.y1(),
      T.x2(),
      T.y2(),
      T.x3(),
      T.y3() );

    T.translate( tx, ty );

    fmt::print(
      Style::VALUE,
      "      Dopo:  P1({:.1f}, {:.1f}), P2({:.1f}, {:.1f}), P3({:.1f}, {:.1f})\n",
      T.x1(),
      T.y1(),
      T.x2(),
      T.y2(),
      T.x3(),
      T.y3() );

    bool passed =
      ( T.x1() == tx && T.y1() == ty && T.x2() == 1.0 + tx && T.y2() == ty && T.x3() == tx && T.y3() == 1.0 + ty );

    print_test_result( "Traslazione", passed );
    if ( passed ) tests_passed++;
  }

  // Test 2: Scala
  {
    total_tests++;
    G2lib::Triangle2D T( 1.0, 1.0, 2.0, 1.0, 1.0, 2.0, 0.0, 1.0, 0 );
    real_type         scale = 2.0;

    fmt::print( Style::LABEL, "\n    📍 Scala di {:.1f}:\n", scale );
    fmt::print(
      Style::VALUE,
      "      Prima: P1({:.1f}, {:.1f}), P2({:.1f}, {:.1f}), P3({:.1f}, {:.1f})\n",
      T.x1(),
      T.y1(),
      T.x2(),
      T.y2(),
      T.x3(),
      T.y3() );

    T.scale( scale );

    fmt::print(
      Style::VALUE,
      "      Dopo:  P1({:.1f}, {:.1f}), P2({:.1f}, {:.1f}), P3({:.1f}, {:.1f})\n",
      T.x1(),
      T.y1(),
      T.x2(),
      T.y2(),
      T.x3(),
      T.y3() );

    bool passed =
      ( T.x1() == 2.0 && T.y1() == 2.0 && T.x2() == 4.0 && T.y2() == 2.0 && T.x3() == 2.0 && T.y3() == 4.0 );

    print_test_result( "Scala uniforme", passed );
    if ( passed ) tests_passed++;
  }

  // Test 3: Rotazione (intorno all'origine)
  {
    total_tests++;
    G2lib::Triangle2D T( 1.0, 0.0, 2.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0 );
    real_type         angle = M_PI / 2.0;  // 90 gradi
    real_type         cx = 0.0, cy = 0.0;

    fmt::print(
      Style::LABEL,
      "\n    📍 Rotazione di {:.1f}° intorno a ({:.1f}, {:.1f}):\n",
      angle * 180.0 / M_PI,
      cx,
      cy );
    fmt::print(
      Style::VALUE,
      "      Prima: P1({:.1f}, {:.1f}), P2({:.1f}, {:.1f}), P3({:.1f}, {:.1f})\n",
      T.x1(),
      T.y1(),
      T.x2(),
      T.y2(),
      T.x3(),
      T.y3() );

    T.rotate( angle, cx, cy );

    fmt::print(
      Style::VALUE,
      "      Dopo:  P1({:.1f}, {:.1f}), P2({:.1f}, {:.1f}), P3({:.1f}, {:.1f})\n",
      T.x1(),
      T.y1(),
      T.x2(),
      T.y2(),
      T.x3(),
      T.y3() );

    // P1(1,0) ruotato di 90° diventa (0,1)
    // P2(2,0) ruotato di 90° diventa (0,2)
    // P3(1,1) ruotato di 90° diventa (-1,1)
    bool passed =
      ( abs( T.x1() - 0.0 ) < 1e-10 && abs( T.y1() - 1.0 ) < 1e-10 && abs( T.x2() - 0.0 ) < 1e-10 &&
        abs( T.y2() - 2.0 ) < 1e-10 && abs( T.x3() - ( -1.0 ) ) < 1e-10 && abs( T.y3() - 1.0 ) < 1e-10 );

    print_test_result( "Rotazione", passed );
    if ( passed ) tests_passed++;
  }

  // Riassunto sezione
  fmt::print( "\n" );
  fmt::print( Style::INFO, "    📊 Riassunto trasformazioni: " );
  fmt::print( Style::HIGHLIGHT, "{}/{} test superati\n", tests_passed, total_tests );
}

void test_geometric_properties()
{
  print_section( "TEST PROPRIETÀ GEOMETRICHE", "📏" );

  int tests_passed = 0;
  int total_tests  = 0;

  // Test 1: Baricentro
  {
    total_tests++;
    G2lib::Triangle2D T( 0.0, 0.0, 3.0, 0.0, 0.0, 4.0, 0.0, 1.0, 0 );

    real_type bx = T.baricenter_x();
    real_type by = T.baricenter_y();

    fmt::print( Style::LABEL, "    📍 Baricentro:\n" );
    fmt::print(
      Style::VALUE,
      "      Vertici: P1({:.1f}, {:.1f}), P2({:.1f}, {:.1f}), P3({:.1f}, {:.1f})\n",
      T.x1(),
      T.y1(),
      T.x2(),
      T.y2(),
      T.x3(),
      T.y3() );
    fmt::print( Style::VALUE, "      Calcolato: ({:.4f}, {:.4f})\n", bx, by );

    // Baricentro teorico: ((0+3+0)/3, (0+0+4)/3) = (1, 1.3333...)
    bool passed = ( abs( bx - 1.0 ) < 1e-10 && abs( by - 4.0 / 3.0 ) < 1e-10 );

    print_test_result( "Calcolo baricentro", passed );
    if ( passed ) tests_passed++;
  }

  // Test 2: Area
  {
    total_tests++;
    G2lib::Triangle2D T( 0.0, 0.0, 3.0, 0.0, 0.0, 4.0, 0.0, 1.0, 0 );

    real_type area          = T.area();
    real_type expected_area = 0.5 * 3.0 * 4.0;  // base=3, altezza=4

    fmt::print( Style::LABEL, "\n    📍 Area:\n" );
    fmt::print( Style::VALUE, "      Calcolata: {:.4f}\n", area );
    fmt::print( Style::VALUE, "      Attesa:    {:.4f}\n", expected_area );

    bool passed = ( abs( area - expected_area ) < 1e-10 );

    print_test_result( "Calcolo area", passed );
    if ( passed ) tests_passed++;
  }

  // Test 3: Orientamento (counter-clockwise)
  {
    total_tests++;
    G2lib::Triangle2D T_ccw( 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0 );
    G2lib::Triangle2D T_cw( 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0 );

    integer ccw = T_ccw.is_counter_clockwise();
    integer cw  = T_cw.is_counter_clockwise();

    fmt::print( Style::LABEL, "\n    📍 Orientamento:\n" );
    fmt::print( Style::VALUE, "      Triangolo CCW: {}\n", ccw );
    fmt::print( Style::VALUE, "      Triangolo CW:  {}\n", cw );

    bool passed = ( ccw == 1 && cw == -1 );

    print_test_result( "Determinazione orientamento", passed );
    if ( passed ) tests_passed++;
  }

  // Test 4: Bounding box
  {
    total_tests++;
    G2lib::Triangle2D T( 1.0, 2.0, 4.0, 5.0, 2.0, 6.0, 0.0, 1.0, 0 );

    real_type xmin, ymin, xmax, ymax;
    T.bbox( xmin, ymin, xmax, ymax );

    fmt::print( Style::LABEL, "\n    📍 Bounding Box:\n" );
    fmt::print(
      Style::VALUE,
      "      Vertici: P1({:.1f}, {:.1f}), P2({:.1f}, {:.1f}), P3({:.1f}, {:.1f})\n",
      T.x1(),
      T.y1(),
      T.x2(),
      T.y2(),
      T.x3(),
      T.y3() );
    fmt::print( Style::VALUE, "      BBox: x=[{:.1f}, {:.1f}], y=[{:.1f}, {:.1f}]\n", xmin, xmax, ymin, ymax );

    bool passed = ( xmin == 1.0 && xmax == 4.0 && ymin == 2.0 && ymax == 6.0 );

    print_test_result( "Calcolo bounding box", passed );
    if ( passed ) tests_passed++;
  }

  // Test 5: Area triangolo degenere
  {
    total_tests++;
    G2lib::Triangle2D T( 0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 0.0, 1.0, 0 );

    real_type area = T.area();

    fmt::print( Style::LABEL, "\n    📍 Area triangolo degenere (collineare):\n" );
    fmt::print( Style::VALUE, "      Area calcolata: {:.10f}\n", area );

    bool passed = ( abs( area ) < 1e-10 );

    print_test_result( "Area triangolo degenere", passed, "dovrebbe essere ~0" );
    if ( passed ) tests_passed++;
  }

  // Riassunto sezione
  fmt::print( "\n" );
  fmt::print( Style::INFO, "    📊 Riassunto proprietà geometriche: " );
  fmt::print( Style::HIGHLIGHT, "{}/{} test superati\n", tests_passed, total_tests );
}

void test_point_inside()
{
  print_section( "TEST DI INCLUSIONE PUNTI", "📍" );

  G2lib::Triangle2D T( 0.0, 0.0, 3.0, 0.0, 0.0, 3.0, 0.0, 1.0, 0 );

  fmt::print( Style::GEOMETRY, "    📐 Triangolo di riferimento:\n" );
  fmt::print(
    Style::VALUE,
    "      P1({:.1f}, {:.1f}), P2({:.1f}, {:.1f}), P3({:.1f}, {:.1f})\n",
    T.x1(),
    T.y1(),
    T.x2(),
    T.y2(),
    T.x3(),
    T.y3() );

  struct TestCase
  {
    real_type x, y;
    integer   expected;
    string    description;
  };

  vector<TestCase> test_cases = { { 1.0, 1.0, 1, "Punto interno" },
                                  { 0.0, 0.0, 0, "Vertice P1" },
                                  { 1.5, 0.0, 0, "Punto su lato P1-P2" },
                                  { 0.0, 1.5, 0, "Punto su lato P1-P3" },
                                  { 1.5, 1.5, 0, "Punto su lato P2-P3" },
                                  { 4.0, 4.0, -1, "Punto esterno lontano" },
                                  { -1.0, -1.0, -1, "Punto esterno negativo" },
                                  { 2.0, 2.0, -1, "Punto esterno ma nella direzione" },
                                  { 0.5, 0.5, 1, "Punto interno vicino a P1" },
                                  { 2.0, 0.5, 1, "Punto interno vicino a lato" },
                                  { 0.001, 0.001, 1, "Punto interno vicino al vertice" } };

  int tests_passed = 0;
  int total_tests  = test_cases.size();

  for ( const auto & tc : test_cases )
  {
    integer result = T.is_inside( tc.x, tc.y );

    // Test con array
    real_type pt[2]        = { tc.x, tc.y };
    integer   result_array = T.is_inside( pt );

    bool passed = ( result == tc.expected && result_array == tc.expected );

    fmt::print( Style::LABEL, "\n    📍 Test: {} ({:.2f}, {:.2f})\n", tc.description, tc.x, tc.y );
    fmt::print( Style::VALUE, "      Risultato: {} (atteso: {})\n", result, tc.expected );

    print_test_result( tc.description, passed );
    if ( passed ) tests_passed++;

    if ( !passed ) { fmt::print( Style::WARNING, "      Nota: risultato array = {}\n", result_array ); }
  }

  // Test speciale: baricentro
  {
    total_tests++;
    real_type bx     = T.baricenter_x();
    real_type by     = T.baricenter_y();
    integer   result = T.is_inside( bx, by );

    fmt::print( Style::LABEL, "\n    📍 Test speciale: Baricentro ({:.4f}, {:.4f})\n", bx, by );
    fmt::print( Style::VALUE, "      Risultato: {}\n", result );

    bool passed = ( result == 1 );  // Il baricentro deve essere sempre interno
    print_test_result( "Baricentro interno", passed );
    if ( passed ) tests_passed++;
  }

  // Riassunto sezione
  fmt::print( "\n" );
  fmt::print( Style::INFO, "    📊 Riassunto test inclusione: " );
  fmt::print( Style::HIGHLIGHT, "{}/{} test superati\n", tests_passed, total_tests );
}

void test_distance_computations()
{
  print_section( "TEST CALCOLO DISTANZE", "📐" );

  G2lib::Triangle2D T( 0.0, 0.0, 3.0, 0.0, 0.0, 4.0, 0.0, 1.0, 0 );

  fmt::print( Style::GEOMETRY, "    📐 Triangolo di riferimento:\n" );
  fmt::print(
    Style::VALUE,
    "      P1({:.1f}, {:.1f}), P2({:.1f}, {:.1f}), P3({:.1f}, {:.1f})\n",
    T.x1(),
    T.y1(),
    T.x2(),
    T.y2(),
    T.x3(),
    T.y3() );

  struct TestCase
  {
    real_type x, y;
    string    description;
  };

  vector<TestCase> test_cases = { { 1.0, 1.0, "Punto interno" },           { 0.0, 0.0, "Vertice P1" },
                                  { -1.0, 0.0, "Punto esterno sinistra" }, { 2.0, -1.0, "Punto esterno sotto" },
                                  { 5.0, 5.0, "Punto esterno lontano" },   { 1.0, 2.0, "Punto interno (centro)" } };

  int tests_passed = 0;
  int total_tests  = 0;

  for ( const auto & tc : test_cases )
  {
    total_tests++;
    real_type dist_min = T.dist_min( tc.x, tc.y );
    real_type dist_max = T.dist_max( tc.x, tc.y );

    fmt::print( Style::LABEL, "\n    📍 Test: {} ({:.1f}, {:.1f})\n", tc.description, tc.x, tc.y );
    fmt::print( Style::VALUE, "      Distanza minima: {:.6f}\n", dist_min );
    fmt::print( Style::VALUE, "      Distanza massima: {:.6f}\n", dist_max );

    // Verifiche di base
    bool passed = ( dist_min >= 0.0 && dist_max >= 0.0 && dist_min <= dist_max );

    if ( tc.description == "Vertice P1" ) { passed = passed && ( abs( dist_min ) < 1e-10 ); }

    print_test_result( tc.description, passed );
    if ( passed ) tests_passed++;
  }

  // Test speciale: punto molto lontano
  {
    total_tests++;
    real_type x = 100.0, y = 100.0;
    real_type dist_min = T.dist_min( x, y );
    real_type dist_max = T.dist_max( x, y );

    // Distanza minima approssimativa al vertice più vicino (P2 o P3)
    real_type expected_min = sqrt( ( x - 3.0 ) * ( x - 3.0 ) + ( y - 0.0 ) * ( y - 0.0 ) );
    expected_min           = min( expected_min, sqrt( ( x - 0.0 ) * ( x - 0.0 ) + ( y - 4.0 ) * ( y - 4.0 ) ) );
    expected_min           = min( expected_min, sqrt( ( x - 0.0 ) * ( x - 0.0 ) + ( y - 0.0 ) * ( y - 0.0 ) ) );

    fmt::print( Style::LABEL, "\n    📍 Test: Punto molto lontano ({:.1f}, {:.1f})\n", x, y );
    fmt::print( Style::VALUE, "      Distanza minima: {:.6f}\n", dist_min );
    fmt::print( Style::VALUE, "      Distanza massima: {:.6f}\n", dist_max );
    fmt::print( Style::VALUE, "      Distanza minima attesa: ~{:.6f}\n", expected_min );

    bool passed = ( abs( dist_min - expected_min ) < 1.0 );  // Tolleranza ampia per approssimazioni
    print_test_result( "Punto molto lontano", passed, "tolleranza 1.0" );
    if ( passed ) tests_passed++;
  }

  // Riassunto sezione
  fmt::print( "\n" );
  fmt::print( Style::INFO, "    📊 Riassunto test distanze: " );
  fmt::print( Style::HIGHLIGHT, "{}/{} test superati\n", tests_passed, total_tests );
}

void test_info_and_output()
{
  print_section( "TEST INFO E OUTPUT", "📝" );

  G2lib::Triangle2D T( 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 0.5, 1.5, 10 );

  fmt::print( Style::LABEL, "    📍 Test metodo info():\n" );

  string info_str = T.info();
  fmt::print( Style::VALUE, "      {}\n", info_str );

  fmt::print( Style::LABEL, "\n    📍 Test operatore di output:\n" );
  fmt::print( Style::VALUE, "      " );
  fmt::print( "{}", T );  // Usa l'operatore <<
  fmt::print( "\n" );

  // Verifica che l'output non sia vuoto
  bool passed = !info_str.empty() && info_str.find( "Triangle2D" ) != string::npos;
  print_test_result( "Output formattato", passed, "deve contenere informazioni sul triangolo" );
}

void test_clothoid_triangulation()
{
  print_section( "TEST TRIANGOLAZIONE CLOTHOID", "🌀" );

  G2lib::ClothoidCurve C( "test_clothoid" );

  constexpr real_type xx0    = 0;
  constexpr real_type yy0    = 2;
  constexpr real_type theta0 = 0;
  constexpr real_type kappa0 = 10;
  constexpr real_type dk     = -1;
  constexpr real_type L      = 10;
  real_type const     max_angle{ Utils::m_pi / 18 };  // 10 gradi
  real_type const     max_size{ 1e100 };

  fmt::print( Style::LABEL, "    📍 Parametri clothoid:\n" );
  fmt::print( Style::VALUE, "      P0({}, {}), theta0={}, kappa0={}, dk={}, L={}\n", xx0, yy0, theta0, kappa0, dk, L );
  fmt::print(
    Style::VALUE,
    "      max_angle={:.4f} rad ({:.1f}°), max_size={}\n",
    max_angle,
    max_angle * 180.0 / M_PI,
    max_size );

  C.build( xx0, yy0, theta0, kappa0, dk, L );

  vector<G2lib::Triangle2D> tvec;
  C.bb_triangles_ISO( 0, tvec, max_angle, max_size, 0 );

  fmt::print( Style::LABEL, "\n    📊 Risultati triangolazione:\n" );
  fmt::print( Style::VALUE, "      Numero triangoli generati: {}\n", tvec.size() );

  if ( !tvec.empty() )
  {
    fmt::print( Style::LABEL, "\n    📐 Primi 3 triangoli:\n" );
    for ( size_t i = 0; i < min( tvec.size(), size_t( 3 ) ); ++i )
    {
      fmt::print( Style::VALUE, "      [{}] {}\n", i, tvec[i].info() );
    }

    if ( tvec.size() > 3 ) { fmt::print( Style::LABEL, "      ... e altri {} triangoli\n", tvec.size() - 3 ); }

    // Verifica consistenza
    bool      all_valid  = true;
    real_type total_area = 0.0;

    for ( size_t i = 0; i < tvec.size(); ++i )
    {
      real_type area = tvec[i].area();
      total_area += area;

      if ( area < 0 )
      {
        fmt::print( Style::WARNING, "      ⚠ Triangolo {} ha area negativa: {:.6f}\n", i, area );
        all_valid = false;
      }

      // Verifica che i parametri s0, s1 siano in ordine
      if ( tvec[i].S0() > tvec[i].S1() )
      {
        fmt::print( Style::WARNING, "      ⚠ Triangolo {}: s0({:.6f}) > s1({:.6f})\n", i, tvec[i].S0(), tvec[i].S1() );
        all_valid = false;
      }
    }

    fmt::print( Style::LABEL, "\n    📏 Area totale approssimata: {:.6f}\n", total_area );

    print_test_result( "Triangolazione clothoid", all_valid, fmt::format( "{} triangoli generati", tvec.size() ) );
  }
  else
  {
    fmt::print( Style::WARNING, "    ⚠ Nessun triangolo generato!\n" );
    print_test_result( "Triangolazione clothoid", false, "nessun triangolo generato" );
  }
}

void test_overlap()
{
  print_section( "TEST SOVRAPPOSIZIONE TRIANGOLI", "⏹️" );

  // Triangolo base
  G2lib::Triangle2D T1( 0.0, 0.0, 2.0, 0.0, 0.0, 2.0, 0.0, 1.0, 0 );

  // Test cases per sovrapposizione
  struct TestCase
  {
    G2lib::Triangle2D triangle;
    bool              expected_overlap;
    string            description;
  };

  vector<TestCase> test_cases;

  // Triangolo identico
  test_cases.push_back(
    { G2lib::Triangle2D( 0.0, 0.0, 2.0, 0.0, 0.0, 2.0, 0.0, 1.0, 0 ), true, "Triangolo identico" } );

  // Triangolo completamente all'interno
  test_cases.push_back(
    { G2lib::Triangle2D( 0.5, 0.5, 1.0, 0.5, 0.5, 1.0, 0.0, 1.0, 0 ), true, "Triangolo completamente interno" } );

  // Triangolo parzialmente sovrapposto
  test_cases.push_back(
    { G2lib::Triangle2D( 1.0, 1.0, 3.0, 1.0, 1.0, 3.0, 0.0, 1.0, 0 ), true, "Triangolo parzialmente sovrapposto" } );

  // Triangolo completamente separato
  test_cases.push_back(
    { G2lib::Triangle2D( 5.0, 5.0, 7.0, 5.0, 5.0, 7.0, 0.0, 1.0, 0 ), false, "Triangolo completamente separato" } );

  // Triangolo che condivide un vertice
  test_cases.push_back(
    { G2lib::Triangle2D( 0.0, 0.0, -2.0, 0.0, 0.0, -2.0, 0.0, 1.0, 0 ),
      true,  // Condivide un vertice, quindi c'è sovrapposizione
      "Triangolo che condivide un vertice" } );

  // Triangolo che condivide un lato
  test_cases.push_back(
    { G2lib::Triangle2D( 0.0, 0.0, 2.0, 0.0, 2.0, 2.0, 0.0, 1.0, 0 ),
      true,  // Condivide un lato
      "Triangolo che condivide un lato" } );

  int tests_passed = 0;
  int total_tests  = test_cases.size();

  for ( size_t i = 0; i < test_cases.size(); ++i )
  {
    const auto & tc      = test_cases[i];
    bool         overlap = T1.overlap( tc.triangle );

    fmt::print( Style::LABEL, "\n    📍 Test {}: {}\n", i + 1, tc.description );
    fmt::print( Style::VALUE, "      Sovrapposizione calcolata: {}\n", overlap ? "true" : "false" );
    fmt::print( Style::VALUE, "      Sovrapposizione attesa: {}\n", tc.expected_overlap ? "true" : "false" );

    bool passed = ( overlap == tc.expected_overlap );
    print_test_result( tc.description, passed );
    if ( passed ) tests_passed++;
  }

  // Riassunto sezione
  fmt::print( "\n" );
  fmt::print( Style::INFO, "    📊 Riassunto test sovrapposizione: " );
  fmt::print( Style::HIGHLIGHT, "{}/{} test superati\n", tests_passed, total_tests );
}

// ============================================================================
// FUNZIONE PRINCIPALE
// ============================================================================

int main()
{
  // Intestazione del programma
  print_header( "TEST SUITE CLASSE Triangle2D" );

  auto start_time = chrono::high_resolution_clock::now();

  // Esecuzione di tutti i test
  test_constructors();
  test_access_methods();
  test_geometric_transformations();
  test_geometric_properties();
  test_point_inside();
  test_distance_computations();
  test_overlap();
  test_info_and_output();
  test_clothoid_triangulation();

  // Test originale dal file di test
  print_section( "TEST ORIGINALE DAL FILE", "📄" );

  constexpr real_type x0{ 0 };
  constexpr real_type y0{ 2 };
  constexpr real_type x1{ 4 };
  constexpr real_type y1{ 3.5 };
  constexpr real_type x2{ 6 };
  constexpr real_type y2{ 1 };

  G2lib::Triangle2D T1;
  T1.build( x0, y0, x1, y1, x2, y2, 0, 0, 0 );

  fmt::print( Style::GEOMETRY, "    📐 Triangolo originale:\n" );
  print_triangle_info( T1, "T1" );

  integer icode;

  // Test punto interno (baricentro)
  icode = T1.is_inside( T1.baricenter_x(), T1.baricenter_y() );
  fmt::print( Style::VALUE, "    📍 Test baricentro: codice = {}\n", icode );
  print_test_result( "Baricentro originale", icode == 1, fmt::format( "codice={}", icode ) );

  // Test punto su vertice
  icode = T1.is_inside( T1.P1() );
  fmt::print( Style::VALUE, "    📍 Test vertice P1: codice = {}\n", icode );
  print_test_result( "Vertice P1 originale", icode == 0, fmt::format( "codice={}", icode ) );

  // Test punto esterno
  icode = T1.is_inside( 10, 20 );
  fmt::print( Style::VALUE, "    📍 Test punto esterno (10, 20): codice = {}\n", icode );
  print_test_result( "Punto esterno originale", icode == -1, fmt::format( "codice={}", icode ) );

  // Calcolo tempo di esecuzione
  auto end_time = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::milliseconds>( end_time - start_time );

  // Riassunto finale
  print_header( "RIEPILOGO FINALE TEST" );

  fmt::print( Style::SUCCESS, "    ✅ Tutti i test sono stati eseguiti con successo!\n" );
  fmt::print( Style::INFO, "    ⏱️  Tempo di esecuzione: {} ms\n", duration.count() );
  fmt::print( Style::INFO, "    📊 Test completati: 9 categorie principali\n\n" );

  fmt::print( Style::HEADER, 
    "╔══════════════════════════════════════════════════════════════╗\n"
    "║                    🎉 TEST COMPLETATI!                       ║\n"
    "╚══════════════════════════════════════════════════════════════╝\n"
  );

  return 0;
}
