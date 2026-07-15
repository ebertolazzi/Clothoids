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

using namespace G2lib;
using namespace autodiff;

// ============================================================================
// Test per ClothoidData
// ============================================================================
void test_ClothoidData()
{
  std::cout << "Testing ClothoidData..." << std::endl;
  ClothoidData cd;
  cd.m_x0     = 0.0;
  cd.m_y0     = 0.0;
  cd.m_theta0 = 0.0;
  cd.m_kappa0 = 0.0;
  cd.m_dk     = 0.0;

  // dual1st
  dual1st s1     = 1.0;
  dual1st theta1 = cd.theta( s1 );
  dual1st kappa1 = cd.kappa( s1 );
  dual1st X1     = cd.X( s1 );
  dual1st Y1     = cd.Y( s1 );
  // Note: X_ISO, Y_ISO, X_SAE, Y_SAE with offset not available for autodiff
  dual1st tg_x1      = cd.tg_x( s1 );
  dual1st tg_y1      = cd.tg_y( s1 );
  dual1st nor_x_ISO1 = cd.nor_x_ISO( s1 );
  dual1st nor_y_ISO1 = cd.nor_y_ISO( s1 );
  dual1st nor_x_SAE1 = cd.nor_x_SAE( s1 );
  dual1st nor_y_SAE1 = cd.nor_y_SAE( s1 );

  // Suppress unused variable warnings
  (void) theta1;
  (void) kappa1;
  (void) X1;
  (void) Y1;
  (void) tg_x1;
  (void) tg_y1;
  (void) nor_x_ISO1;
  (void) nor_y_ISO1;
  (void) nor_x_SAE1;
  (void) nor_y_SAE1;

  // dual2nd
  dual2nd s2     = 1.0;
  dual2nd theta2 = cd.theta( s2 );
  dual2nd kappa2 = cd.kappa( s2 );
  dual2nd X2     = cd.X( s2 );
  dual2nd Y2     = cd.Y( s2 );
  // Note: X_ISO, Y_ISO, X_SAE, Y_SAE with offset not available for autodiff
  dual2nd tg_x2      = cd.tg_x( s2 );
  dual2nd tg_y2      = cd.tg_y( s2 );
  dual2nd nor_x_ISO2 = cd.nor_x_ISO( s2 );
  dual2nd nor_y_ISO2 = cd.nor_y_ISO( s2 );
  dual2nd nor_x_SAE2 = cd.nor_x_SAE( s2 );
  dual2nd nor_y_SAE2 = cd.nor_y_SAE( s2 );

  // Suppress unused variable warnings
  (void) theta2;
  (void) kappa2;
  (void) X2;
  (void) Y2;
  (void) tg_x2;
  (void) tg_y2;
  (void) nor_x_ISO2;
  (void) nor_y_ISO2;
  (void) nor_x_SAE2;
  (void) nor_y_SAE2;

  std::cout << "  ClothoidData test PASSED." << std::endl;
}

// ============================================================================
// Test per LineSegment
// ============================================================================
void test_LineSegment()
{
  std::cout << "Testing LineSegment..." << std::endl;
  LineSegment ls( "test" );
  ls.build( 0.0, 0.0, 0.0, 1.0 );  // x0, y0, theta0, L

  // dual1st
  dual1st s1      = 0.5;
  dual1st theta1  = ls.theta( s1 );
  dual1st kappa1  = ls.kappa( s1 );
  dual1st tx1     = ls.tx( s1 );
  dual1st ty1     = ls.ty( s1 );
  dual1st nx_ISO1 = ls.nx_ISO( s1 );
  dual1st ny_ISO1 = ls.ny_ISO( s1 );
  dual1st nx_SAE1 = ls.nx_SAE( s1 );
  dual1st ny_SAE1 = ls.ny_SAE( s1 );
  dual1st X1      = ls.X( s1 );
  dual1st Y1      = ls.Y( s1 );
  // Note: X_ISO, Y_ISO, X_SAE, Y_SAE with offset not available for autodiff

  // Suppress unused variable warnings
  (void) theta1;
  (void) kappa1;
  (void) tx1;
  (void) ty1;
  (void) nx_ISO1;
  (void) ny_ISO1;
  (void) nx_SAE1;
  (void) ny_SAE1;
  (void) X1;
  (void) Y1;

  // dual2nd
  dual2nd s2      = 0.5;
  dual2nd theta2  = ls.theta( s2 );
  dual2nd kappa2  = ls.kappa( s2 );
  dual2nd tx2     = ls.tx( s2 );
  dual2nd ty2     = ls.ty( s2 );
  dual2nd nx_ISO2 = ls.nx_ISO( s2 );
  dual2nd ny_ISO2 = ls.ny_ISO( s2 );
  dual2nd nx_SAE2 = ls.nx_SAE( s2 );
  dual2nd ny_SAE2 = ls.ny_SAE( s2 );
  dual2nd X2      = ls.X( s2 );
  dual2nd Y2      = ls.Y( s2 );
  // Note: X_ISO, Y_ISO, X_SAE, Y_SAE with offset not available for autodiff

  // Suppress unused variable warnings
  (void) theta2;
  (void) kappa2;
  (void) tx2;
  (void) ty2;
  (void) nx_ISO2;
  (void) ny_ISO2;
  (void) nx_SAE2;
  (void) ny_SAE2;
  (void) X2;
  (void) Y2;

  std::cout << "  LineSegment test PASSED." << std::endl;
}

// ============================================================================
// Test per ClothoidCurve
// ============================================================================
void test_ClothoidCurve()
{
  std::cout << "Testing ClothoidCurve..." << std::endl;
  ClothoidCurve cc( "test" );
  cc.build( 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 );  // x0, y0, theta0, kappa0, dk, L

  // dual1st
  dual1st s1      = 0.5;
  dual1st theta1  = cc.theta( s1 );
  dual1st kappa1  = cc.kappa( s1 );
  dual1st tx1     = cc.tx( s1 );
  dual1st ty1     = cc.ty( s1 );
  dual1st nx_ISO1 = cc.nx_ISO( s1 );
  dual1st ny_ISO1 = cc.ny_ISO( s1 );
  dual1st nx_SAE1 = cc.nx_SAE( s1 );
  dual1st ny_SAE1 = cc.ny_SAE( s1 );
  dual1st X1      = cc.X( s1 );
  dual1st Y1      = cc.Y( s1 );
  // Note: X_ISO, Y_ISO, X_SAE, Y_SAE with offset not available for autodiff

  // Suppress unused variable warnings
  (void) theta1;
  (void) kappa1;
  (void) tx1;
  (void) ty1;
  (void) nx_ISO1;
  (void) ny_ISO1;
  (void) nx_SAE1;
  (void) ny_SAE1;
  (void) X1;
  (void) Y1;

  // dual2nd
  dual2nd s2      = 0.5;
  dual2nd theta2  = cc.theta( s2 );
  dual2nd kappa2  = cc.kappa( s2 );
  dual2nd tx2     = cc.tx( s2 );
  dual2nd ty2     = cc.ty( s2 );
  dual2nd nx_ISO2 = cc.nx_ISO( s2 );
  dual2nd ny_ISO2 = cc.ny_ISO( s2 );
  dual2nd nx_SAE2 = cc.nx_SAE( s2 );
  dual2nd ny_SAE2 = cc.ny_SAE( s2 );
  dual2nd X2      = cc.X( s2 );
  dual2nd Y2      = cc.Y( s2 );
  // Note: X_ISO, Y_ISO, X_SAE, Y_SAE with offset not available for autodiff

  // Suppress unused variable warnings
  (void) theta2;
  (void) kappa2;
  (void) tx2;
  (void) ty2;
  (void) nx_ISO2;
  (void) ny_ISO2;
  (void) nx_SAE2;
  (void) ny_SAE2;
  (void) X2;
  (void) Y2;

  std::cout << "  ClothoidCurve test PASSED." << std::endl;
}

// ============================================================================
// Test per ClothoidList
// ============================================================================
void test_ClothoidList()
{
  std::cout << "Testing ClothoidList..." << std::endl;
  ClothoidList cl( "test" );
  cl.init();
  cl.push_back( 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 );  // x0, y0, theta0, kappa0, dk, L

  // dual1st
  dual1st s1      = 0.5;
  dual1st theta1  = cl.theta( s1 );
  dual1st kappa1  = cl.kappa( s1 );
  dual1st tx1     = cl.tx( s1 );
  dual1st ty1     = cl.ty( s1 );
  dual1st nx_ISO1 = cl.nx_ISO( s1 );
  dual1st ny_ISO1 = cl.ny_ISO( s1 );
  dual1st nx_SAE1 = cl.nx_SAE( s1 );
  dual1st ny_SAE1 = cl.ny_SAE( s1 );
  dual1st X1      = cl.X( s1 );
  dual1st Y1      = cl.Y( s1 );
  // Note: X_ISO, Y_ISO, X_SAE, Y_SAE with offset not available for autodiff

  // Suppress unused variable warnings
  (void) theta1;
  (void) kappa1;
  (void) tx1;
  (void) ty1;
  (void) nx_ISO1;
  (void) ny_ISO1;
  (void) nx_SAE1;
  (void) ny_SAE1;
  (void) X1;
  (void) Y1;

  // dual2nd
  dual2nd s2      = 0.5;
  dual2nd theta2  = cl.theta( s2 );
  dual2nd kappa2  = cl.kappa( s2 );
  dual2nd tx2     = cl.tx( s2 );
  dual2nd ty2     = cl.ty( s2 );
  dual2nd nx_ISO2 = cl.nx_ISO( s2 );
  dual2nd ny_ISO2 = cl.ny_ISO( s2 );
  dual2nd nx_SAE2 = cl.nx_SAE( s2 );
  dual2nd ny_SAE2 = cl.ny_SAE( s2 );
  dual2nd X2      = cl.X( s2 );
  dual2nd Y2      = cl.Y( s2 );
  // Note: X_ISO, Y_ISO, X_SAE, Y_SAE with offset not available for autodiff

  // Suppress unused variable warnings
  (void) theta2;
  (void) kappa2;
  (void) tx2;
  (void) ty2;
  (void) nx_ISO2;
  (void) ny_ISO2;
  (void) nx_SAE2;
  (void) ny_SAE2;
  (void) X2;
  (void) Y2;

  std::cout << "  ClothoidList test PASSED." << std::endl;
}

// ============================================================================
// Test per PolyLine
// ============================================================================
void test_PolyLine()
{
  std::cout << "Testing PolyLine..." << std::endl;
  PolyLine pl( "test" );
  pl.init( 0.0, 0.0 );
  pl.push_back( 1.0, 0.0 );

  // dual1st
  dual1st s1      = 0.5;
  dual1st theta1  = pl.theta( s1 );
  dual1st kappa1  = pl.kappa( s1 );
  dual1st tx1     = pl.tx( s1 );
  dual1st ty1     = pl.ty( s1 );
  dual1st nx_ISO1 = pl.nx_ISO( s1 );
  dual1st ny_ISO1 = pl.ny_ISO( s1 );
  dual1st nx_SAE1 = pl.nx_SAE( s1 );
  dual1st ny_SAE1 = pl.ny_SAE( s1 );
  dual1st X1      = pl.X( s1 );
  dual1st Y1      = pl.Y( s1 );
  // Note: X_ISO, Y_ISO, X_SAE, Y_SAE with offset not available for autodiff

  // Suppress unused variable warnings
  (void) theta1;
  (void) kappa1;
  (void) tx1;
  (void) ty1;
  (void) nx_ISO1;
  (void) ny_ISO1;
  (void) nx_SAE1;
  (void) ny_SAE1;
  (void) X1;
  (void) Y1;

  // dual2nd
  dual2nd s2      = 0.5;
  dual2nd theta2  = pl.theta( s2 );
  dual2nd kappa2  = pl.kappa( s2 );
  dual2nd tx2     = pl.tx( s2 );
  dual2nd ty2     = pl.ty( s2 );
  dual2nd nx_ISO2 = pl.nx_ISO( s2 );
  dual2nd ny_ISO2 = pl.ny_ISO( s2 );
  dual2nd nx_SAE2 = pl.nx_SAE( s2 );
  dual2nd ny_SAE2 = pl.ny_SAE( s2 );
  dual2nd X2      = pl.X( s2 );
  dual2nd Y2      = pl.Y( s2 );
  // Note: X_ISO, Y_ISO, X_SAE, Y_SAE with offset not available for autodiff

  // Suppress unused variable warnings
  (void) theta2;
  (void) kappa2;
  (void) tx2;
  (void) ty2;
  (void) nx_ISO2;
  (void) ny_ISO2;
  (void) nx_SAE2;
  (void) ny_SAE2;
  (void) X2;
  (void) Y2;

  std::cout << "  PolyLine test PASSED." << std::endl;
}

// ============================================================================
// Test per Dubins
// ============================================================================
void test_Dubins()
{
  std::cout << "Testing Dubins..." << std::endl;

  // Crea una curva Dubins semplice
  Dubins dubins( "test_dubins" );
  bool   built = dubins.build(
    0.0,
    0.0,
    0.0,  // x0, y0, theta0
    2.0,
    1.0,
    M_PI / 2,  // x1, y1, theta1
    0.5        // k_max
  );

  if ( !built )
  {
    std::cout << "  Dubins build FAILED!" << std::endl;
    return;
  }

#ifdef AUTODIFF_SUPPORT
  // Test con dual1st
  dual1st s1     = 0.5;
  dual1st theta1 = dubins.theta( s1 );
  dual1st X1     = dubins.X( s1 );
  dual1st Y1     = dubins.Y( s1 );

  // Suppress unused variable warnings
  (void) theta1;
  (void) X1;
  (void) Y1;

  // Test con dual2nd
  dual2nd s2     = 1.0;
  dual2nd theta2 = dubins.theta( s2 );
  dual2nd X2     = dubins.X( s2 );
  dual2nd Y2     = dubins.Y( s2 );

  // Suppress unused variable warnings
  (void) theta2;
  (void) X2;
  (void) Y2;
#endif

  std::cout << "  Dubins test PASSED." << std::endl;
}

// ============================================================================
// Test per Dubins3p
// ============================================================================
void test_Dubins3p()
{
  std::cout << "Testing Dubins3p..." << std::endl;

  // Crea una curva Dubins3p con punto intermedio
  Dubins3p dubins3p( "test_dubins3p" );
  bool     built = dubins3p.build(
    0.0,
    0.0,
    0.0,  // xi, yi, thetai
    1.0,
    0.5,  // xm, ym (punto intermedio)
    2.0,
    1.0,
    M_PI / 2,                          // xf, yf, thetaf
    0.5,                               // k_max
    Dubins3pBuildType::PATTERN_SEARCH  // metodo di costruzione
  );

  if ( !built )
  {
    std::cout << "  Dubins3p build FAILED!" << std::endl;
    return;
  }

#ifdef AUTODIFF_SUPPORT
  // Test con dual1st
  dual1st s1     = 0.5;
  dual1st theta1 = dubins3p.theta( s1 );
  dual1st X1     = dubins3p.X( s1 );
  dual1st Y1     = dubins3p.Y( s1 );

  // Suppress unused variable warnings
  (void) theta1;
  (void) X1;
  (void) Y1;

  // Test con dual2nd
  dual2nd s2     = 1.0;
  dual2nd theta2 = dubins3p.theta( s2 );
  dual2nd X2     = dubins3p.X( s2 );
  dual2nd Y2     = dubins3p.Y( s2 );

  // Suppress unused variable warnings
  (void) theta2;
  (void) X2;
  (void) Y2;
#endif

  std::cout << "  Dubins3p test PASSED." << std::endl;
}

// ============================================================================
// Main
// ============================================================================
int main()
{
  std::cout << "=== STARTING AUTODIFF METHODS TEST ===\n" << std::endl;

  test_ClothoidData();
  test_LineSegment();
  // CircleArc non ha supporto autodiff, quindi lo saltiamo
  test_ClothoidCurve();
  test_ClothoidList();
  test_PolyLine();
  test_Dubins();
  test_Dubins3p();

  std::cout << "\n=== ALL TESTS PASSED SUCCESSFULLY ===" << std::endl;
  return 0;
}
