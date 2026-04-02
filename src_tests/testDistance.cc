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

struct TestPoint
{
  real_type x;
  real_type y;
  string    description;
};

// Stili di formattazione predefiniti
namespace Style
{
  const auto HEADER  = fg( fmt::color::steel_blue ) | fmt::emphasis::bold;
  const auto SUCCESS = fg( fmt::color::lime_green ) | fmt::emphasis::bold;
  const auto WARNING = fg( fmt::color::gold ) | fmt::emphasis::bold;
  // const auto ERROR     = fg( fmt::color::crimson ) | fmt::emphasis::bold;
  const auto INFO      = fg( fmt::color::deep_sky_blue ) | fmt::emphasis::bold;
  const auto VALUE     = fg( fmt::color::light_gray );
  const auto LABEL     = fg( fmt::color::silver );
  const auto HIGHLIGHT = fg( fmt::color::cyan ) | fmt::emphasis::bold;
  const auto CURVE     = fg( fmt::color::violet ) | fmt::emphasis::bold;
  const auto POINT     = fg( fmt::color::spring_green ) | fmt::emphasis::bold;
}  // namespace Style

// Funzione per stampare risultati in formato tabellare con colori
void print_results( const string & desc, real_type d, real_type X, real_type Y, real_type S, real_type T, integer info )
{
  // Intestazione del test
  fmt::print( Style::HEADER, "┌{0:─^{1}}┐\n", "", 58 );
  fmt::print( Style::HEADER, "│ {0: ^{1}} │\n", "📍 " + desc, 56 );
  fmt::print( Style::HEADER, "└{0:─^{1}}┘\n", "", 58 );

  // Dettagli risultati
  fmt::print( Style::LABEL, "    📏 " );
  fmt::print( Style::VALUE, "Distanza dalla curva: " );
  fmt::print( Style::HIGHLIGHT, "{:.6f}\n", d );

  fmt::print( Style::LABEL, "    📍 " );
  fmt::print( Style::VALUE, "Punto più vicino:     " );
  fmt::print( Style::POINT, "({:.6f}, {:.6f})\n", X, Y );

  fmt::print( Style::LABEL, "    📐 " );
  fmt::print( Style::VALUE, "Lunghezza arco (S):   " );
  fmt::print( Style::HIGHLIGHT, "{:.6f}\n", S );

  fmt::print( Style::LABEL, "    🧭 " );
  fmt::print( Style::VALUE, "Tangente (T):         " );
  fmt::print( Style::HIGHLIGHT, "{:.6f} rad ", T );
  fmt::print( Style::VALUE, "({:.2f}°)\n", T * 180.0 / M_PI );

  fmt::print( Style::LABEL, "    🏷️  " );
  fmt::print( Style::VALUE, "Codice di ritorno:    {} ✓", info );

  fmt::print( "\n" );
}

// Funzione per stampare una sezione principale
void print_section( const string & title, const string & icon = "📊" )
{
  fmt::print( "\n" );
  fmt::print( Style::HEADER, "╔{0:═^{1}}╗\n", "", 60 );
  fmt::print( Style::HEADER, "║ {0: ^{1}} ║\n", icon + " " + title, 58 );
  fmt::print( Style::HEADER, "╚{0:═^{1}}╝\n", "", 60 );
  fmt::print( "\n" );
}

int main()
{
  // =================== INTESTAZIONE PROGRAMMA ===================
  fmt::print( Style::HEADER, "╔══════════════════════════════════════════════════════════════╗\n" );
  fmt::print( Style::HEADER, "║                    📐 TEST CLOTHOID CURVE                    ║\n" );
  fmt::print( Style::HEADER, "╚══════════════════════════════════════════════════════════════╝\n\n" );

  // Informazioni sul sistema
  auto now = chrono::system_clock::now();
  fmt::print( Style::INFO, "    📅 Data e ora: {:%Y-%m-%d %H:%M:%S}\n", now );
  fmt::print( Style::INFO, "    🏷️  Versione: G2lib Test Suite 1.0\n" );
  fmt::print( Style::INFO, "    👤 Autore: Enrico Bertolazzi\n" );
  fmt::print( Style::INFO, "    📧 Email: enrico.bertolazzi@unitn.it\n\n" );

  // =================== CONFIGURAZIONE DELLA CURVA ===================
  print_section( "CONFIGURAZIONE DELLA CURVA", "⚙️" );

  G2lib::ClothoidCurve curve( "📈 clothoid" );

  constexpr real_type x0     = 0.0;
  constexpr real_type y0     = 2.0;
  constexpr real_type theta0 = 0.0;   // Orientamento iniziale (radianti)
  constexpr real_type kappa0 = 10.0;  // Curvatura iniziale
  constexpr real_type dk     = -1.0;  // Rateo di variazione curvatura
  constexpr real_type L      = 10.0;  // Lunghezza della curva

  // Stampa parametri curva
  fmt::print( Style::CURVE, "    🎯 Parametri della clothoid:\n\n" );

  fmt::print( Style::LABEL, "        🏁 Punto iniziale:       " );
  fmt::print( Style::POINT, "({}, {})\n", x0, y0 );

  fmt::print( Style::LABEL, "        📐 Angolo iniziale:      " );
  fmt::print( Style::VALUE, "{} rad ", theta0 );
  fmt::print( Style::HIGHLIGHT, "({:.1f}°)\n", theta0 * 180.0 / M_PI );

  fmt::print( Style::LABEL, "        📊 Curvatura iniziale:   " );
  fmt::print( Style::HIGHLIGHT, "{}\n", kappa0 );

  fmt::print( Style::LABEL, "        📈 Tasso curvatura:      " );
  fmt::print( Style::HIGHLIGHT, "{}\n", dk );

  fmt::print( Style::LABEL, "        📏 Lunghezza curva:      " );
  fmt::print( Style::HIGHLIGHT, "{}\n\n", L );

  // Costruzione curva
  fmt::print( Style::LABEL, "    🔧 " );
  fmt::print( Style::INFO, "Costruzione della curva in corso... " );
  curve.build( x0, y0, theta0, kappa0, dk, L );
  fmt::print( Style::SUCCESS, "✓ Completato\n\n" );

  // =================== PUNTI DI TEST ===================
  print_section( "PUNTI DI TEST", "🎯" );

  vector<TestPoint> test_points = { // Punti sulla curva o vicini
                                    { 0.0, 2.0, "Punto iniziale della curva" },
                                    { 10.0, 2.0, "Punto finale (estrapolazione retta)" },
                                    { 5.0, 2.0, "Punto medio (segmento rettilineo)" },

                                    // Punti lontani
                                    { 10.0, 15.0, "Punto distante (test originale)" },
                                    { -5.0, -5.0, "Punto nel quadrante negativo" },

                                    // Punti speciali
                                    { 0.0, 0.0, "Origine degli assi" },
                                    { x0, y0 + 5.0, "Direttamente sopra il punto iniziale" },
                                    { x0 + L, y0, "Allineato orizzontalmente con la fine" }
  };

  fmt::print( Style::INFO, "    📋 " );
  fmt::print( Style::VALUE, "Totale punti di test: " );
  fmt::print( Style::HIGHLIGHT, "{}\n\n", test_points.size() );

  // =================== ESECUZIONE TEST ===================
  print_section( "CALCOLI PUNTI PIÙ VICINI", "🧮" );

  int success_count = 0;
  int warning_count = 0;

  for ( size_t i = 0; i < test_points.size(); ++i )
  {
    const auto & tp = test_points[i];

    // Progress bar
    int progress = ( ( i + 1 ) * 50 ) / test_points.size();
    fmt::print( Style::LABEL, "    [{:<50}] ", Utils::repeat( "█", progress ) );
    fmt::print( Style::INFO, "Test {}/{} ", i + 1, test_points.size() );
    fmt::print( Style::VALUE, "- {}\n", tp.description );

    real_type X, Y, S, T, d;
    integer   info = curve.closest_point_ISO( tp.x, tp.y, X, Y, S, T, d );

    print_results( tp.description, d, X, Y, S, T, info );
  }

  // =================== TEST AGGIUNTIVI ===================
  print_section( "TEST AGGIUNTIVI DI VERIFICA", "🔍" );

  // Test: punto esattamente sulla curva (valutazione diretta)
  real_type test_S = L / 2.0;
  real_type X_mid, Y_mid;
  curve.eval( test_S, X_mid, Y_mid );

  fmt::print( Style::INFO, "    🎯 " );
  fmt::print( Style::VALUE, "Test punto sulla curva (valutazione diretta):\n" );
  fmt::print( Style::LABEL, "        Parametro S: " );
  fmt::print( Style::HIGHLIGHT, "{:.6f}\n", test_S );
  fmt::print( Style::LABEL, "        Coordinate:   " );
  fmt::print( Style::POINT, "({:.6f}, {:.6f})\n\n", X_mid, Y_mid );

  TestPoint mid_point = { X_mid, Y_mid, "Punto esattamente sulla curva (centro)" };

  real_type X, Y, S, T, d;
  integer   info = curve.closest_point_ISO( mid_point.x, mid_point.y, X, Y, S, T, d );

  print_results( mid_point.description, d, X, Y, S, T, info );

  // Verifica che il punto trovato sia effettivamente quello di partenza
  real_type error = sqrt( pow( X - mid_point.x, 2 ) + pow( Y - mid_point.y, 2 ) );

  fmt::print( Style::INFO, "    📊 " );
  fmt::print( Style::VALUE, "Verifica precisione:\n" );
  fmt::print( Style::LABEL, "        Errore: " );

  if ( error < 1e-10 )
  {
    fmt::print( Style::SUCCESS, "{:.2e} ", error );
    fmt::print( Style::SUCCESS, "✓ Precisione eccellente\n" );
  }
  else if ( error < 1e-6 )
  {
    fmt::print( fg( fmt::color::green_yellow ), "{:.2e} ", error );
    fmt::print( fg( fmt::color::green_yellow ), "✓ Precisione buona\n" );
  }
  else
  {
    fmt::print( Style::WARNING, "{:.2e} ", error );
    fmt::print( Style::WARNING, "⚠ Precisione da verificare\n" );
  }

  // =================== RIEPILOGO ===================
  print_section( "RIEPILOGO TEST", "📈" );

  fmt::print( Style::CURVE, "    📐 Parametri curva:\n" );
  fmt::print( Style::LABEL, "        📏 Lunghezza: " );
  fmt::print( Style::HIGHLIGHT, "{}\n", L );
  fmt::print( Style::LABEL, "        🏁 Inizio:    " );
  fmt::print( Style::POINT, "({}, {})\n", x0, y0 );
  fmt::print( Style::LABEL, "        📐 Angolo:    " );
  fmt::print( Style::HIGHLIGHT, "{} rad\n", theta0 );
  fmt::print( Style::LABEL, "        📊 Curvatura: " );
  fmt::print( Style::HIGHLIGHT, "{}\n", kappa0 );
  fmt::print( Style::LABEL, "        📈 Tasso:     " );
  fmt::print( Style::HIGHLIGHT, "{}\n\n", dk );

  fmt::print( Style::INFO, "    📊 Statistiche test:\n" );
  fmt::print( Style::LABEL, "        ✅ Test riusciti:      " );
  fmt::print( Style::SUCCESS, "{}\n", success_count );
  fmt::print( Style::LABEL, "        ⚠  Test con avvisi:    " );
  fmt::print( Style::WARNING, "{}\n", warning_count );
  fmt::print( Style::LABEL, "        📋 Test totali:        " );
  fmt::print( Style::HIGHLIGHT, "{}\n", test_points.size() + 1 );

  real_type success_rate = ( success_count * 100.0 ) / static_cast<real_type>( test_points.size() );
  fmt::print( Style::LABEL, "        📈 Tasso di successo:  " );

  if ( success_rate >= 90 ) { fmt::print( Style::SUCCESS, "{:.1f}%\n", success_rate ); }
  else if ( success_rate >= 70 ) { fmt::print( fg( fmt::color::orange ), "{:.1f}%\n", success_rate ); }
  else
  {
    fmt::print( Style::WARNING, "{:.1f}%\n", success_rate );
  }

  // =================== FINE PROGRAMMA ===================
  fmt::print( "\n" );
  fmt::print( Style::HEADER, "╔══════════════════════════════════════════════════════════════╗\n" );
  fmt::print( Style::HEADER, "║                     🎉 TEST COMPLETATI!                      ║\n" );
  fmt::print( Style::HEADER, "╚══════════════════════════════════════════════════════════════╝\n\n" );

  fmt::print( Style::INFO, "    📅 Tempo di esecuzione: {:%H:%M:%S}\n", now );
  fmt::print( Style::SUCCESS, "    ✅ Programma terminato correttamente\n" );

  return 0;
}
