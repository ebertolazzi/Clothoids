/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                      |
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
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include "Fresnel.hh"

#include <cmath>
#include <cfloat>

#define A_THRESOLD   0.01
#define A_SERIE_SIZE 3

namespace G2lib {

  using namespace std ;

  /*
  // This function calculates the fresnel cosine and sine integrals.
  // Input:
  // y = values for which fresnel integrals have to be evaluated  
  //
  // Output:
  // FresnelC = fresnel cosine integral of y
  // FresnelS = fresnel sine integral of y  
  //
  // Adapted from:
  // Atlas for computing mathematical functions : an illustrated guide for
  // practitioners, with programs in C and Mathematica / William J. Thompson.
  // New York : Wiley, c1997.
  //
  // Author: Venkata Sivakanth Telasula
  // email: sivakanth.telasula@gmail.com
  // date: August 11, 2005
  */
  //! \cond NODOC
  static const valueType fn[] = { 0.49999988085884732562,
                                  1.3511177791210715095,
                                  1.3175407836168659241,
                                  1.1861149300293854992,
                                  0.7709627298888346769,
                                  0.4173874338787963957,
                                  0.19044202705272903923,
                                  0.06655998896627697537,
                                  0.022789258616785717418,
                                  0.0040116689358507943804,
                                  0.0012192036851249883877 } ;

  static const valueType fd[] = { 1.0,
                                  2.7022305772400260215,
                                  4.2059268151438492767,
                                  4.5221882840107715516,
                                  3.7240352281630359588,
                                  2.4589286254678152943,
                                  1.3125491629443702962,
                                  0.5997685720120932908,
                                  0.20907680750378849485,
                                  0.07159621634657901433,
                                  0.012602969513793714191,
                                  0.0038302423512931250065 } ;

  static const valueType gn[] = { 0.50000014392706344801,
                                  0.032346434925349128728,
                                  0.17619325157863254363,
                                  0.038606273170706486252,
                                  0.023693692309257725361,
                                  0.007092018516845033662,
                                  0.0012492123212412087428,
                                  0.00044023040894778468486,
                                 -8.80266827476172521e-6,
                                 -1.4033554916580018648e-8,
                                  2.3509221782155474353e-10 } ;

  static const valueType gd[] = { 1.0,
                                  2.0646987497019598937,
                                  2.9109311766948031235,
                                  2.6561936751333032911,
                                  2.0195563983177268073,
                                  1.1167891129189363902,
                                  0.57267874755973172715,
                                  0.19408481169593070798,
                                  0.07634808341431248904,
                                  0.011573247407207865977,
                                  0.0044099273693067311209,
                                 -0.00009070958410429993314 } ;

  static const valueType m_pi        = 3.14159265358979323846264338328  ; // pi
  static const valueType m_pi_2      = 1.57079632679489661923132169164  ; // pi/2
  static const valueType m_1_pi      = 0.318309886183790671537767526745 ; // 1/pi
  static const valueType m_1_sqrt_pi = 0.564189583547756286948079451561 ; // 1/sqrt(pi)

  //! \endcond

  /*
  //  #######                                           
  //  #       #####  ######  ####  #    # ###### #      
  //  #       #    # #      #      ##   # #      #      
  //  #####   #    # #####   ####  # #  # #####  #      
  //  #       #####  #           # #  # # #      #      
  //  #       #   #  #      #    # #   ## #      #      
  //  #       #    # ######  ####  #    # ###### ###### 
  */

  void
  FresnelCS( valueType y, valueType & C, valueType & S ) {
    /*=======================================================*\
      Purpose: This program computes the Fresnel integrals 
               C(x) and S(x) using subroutine FCS
      Input :  x --- Argument of C(x) and S(x)
      Output:  C --- C(x)
               S --- S(x)
      Example:
              x          C(x)          S(x)
             -----------------------------------
             0.0      .00000000      .00000000
             0.5      .49234423      .06473243
             1.0      .77989340      .43825915
             1.5      .44526118      .69750496
             2.0      .48825341      .34341568
             2.5      .45741301      .61918176

      Purpose: Compute Fresnel integrals C(x) and S(x)
      Input :  x --- Argument of C(x) and S(x)
      Output:  C --- C(x)
               S --- S(x)
    \*=======================================================*/

    valueType const eps = 1E-15 ;    
    valueType const x   = y > 0 ? y : -y ;

    if ( x < 1.0 ) {
      valueType twofn, fact, denterm, numterm, sum, term ;

      valueType const s = m_pi_2*(x*x) ;
      valueType const t = -s*s ;

      // Cosine integral series
      twofn   =  0.0 ;
      fact    =  1.0 ;
      denterm =  1.0 ;
      numterm =  1.0 ;
      sum     =  1.0 ;
      do {
        twofn   += 2.0 ;
        fact    *= twofn*(twofn-1.0);
        denterm += 4.0 ;
        numterm *= t ;
        term     = numterm/(fact*denterm) ;
        sum     += term ;
      } while ( std::abs(term) > eps*std::abs(sum) ) ;

      C = x*sum;

      // Sine integral series
      twofn   = 1.0 ;
      fact    = 1.0 ;
      denterm = 3.0 ;
      numterm = 1.0 ;
      sum     = 1.0/3.0 ;
      do {
        twofn   += 2.0 ;
        fact    *= twofn*(twofn-1.0) ;
        denterm += 4.0 ;
        numterm *= t ;
        term     = numterm/(fact*denterm) ;
        sum     += term ;
      } while ( std::abs(term) > eps*std::abs(sum) ) ;

      S = m_pi_2*sum*(x*x*x) ;

    } else if ( x < 6.0 ) {

      // Rational approximation for f
      valueType sumn = 0.0 ;
      valueType sumd = fd[11] ;
      for ( indexType k=10 ; k >= 0 ; --k ) {
        sumn = fn[k] + x*sumn ;
        sumd = fd[k] + x*sumd ;
      }
      valueType f = sumn/sumd ;

      // Rational approximation for g
      sumn = 0.0 ;
      sumd = gd[11] ;
      for ( indexType k=10 ; k >= 0 ; --k ) {
        sumn = gn[k] + x*sumn ;
        sumd = gd[k] + x*sumd ;
      }
      valueType g = sumn/sumd ;

      valueType U    = m_pi_2*(x*x) ;
      valueType SinU = sin(U) ;
      valueType CosU = cos(U) ;
      C = 0.5 + f*SinU - g*CosU ;
      S = 0.5 - f*CosU - g*SinU ;

    } else {

      valueType absterm ;

      // x >= 6; asymptotic expansions for  f  and  g

      valueType const s = m_pi*x*x ;
      valueType const t = -1/(s*s) ;

      // Expansion for f
      valueType numterm = -1.0 ;
      valueType term    =  1.0 ;
      valueType sum     =  1.0 ;
      valueType oldterm =  1.0 ;
      valueType eps10   =  0.1 * eps ;

      do {
        numterm += 4.0 ;
        term    *= numterm*(numterm-2.0)*t ;
        sum     += term ;
        absterm  = std::abs(term) ;
        G2LIB_ASSERT( oldterm >= absterm,
                      "In FresnelCS f not converged to eps, x = " << x <<
                      " oldterm = " << oldterm << " absterm = " << absterm ) ;
        oldterm  = absterm ;
      } while ( absterm > eps10 * std::abs(sum) ) ;

      valueType f = sum / (m_pi*x) ;

      //  Expansion for  g
      numterm = -1.0 ;
      term    =  1.0 ;
      sum     =  1.0 ;
      oldterm =  1.0 ;

      do {
        numterm += 4.0 ;
        term    *= numterm*(numterm+2.0)*t ;
        sum     += term ;
        absterm  = std::abs(term) ;
        G2LIB_ASSERT( oldterm >= absterm,
                      "In FresnelCS g not converged to eps, x = " << x <<
                      " oldterm = " << oldterm << " absterm = " << absterm ) ;
        oldterm  = absterm ;
      } while ( absterm > eps10 * std::abs(sum) ) ;

      valueType g = m_pi*x ; g = sum/(g*g*x) ;

      valueType U    = m_pi_2*(x*x) ;
      valueType SinU = sin(U) ;
      valueType CosU = cos(U) ;
      C = 0.5 + f*SinU - g*CosU ;
      S = 0.5 - f*CosU - g*SinU ;
      
    }
    if ( y < 0 ) { C = -C ; S = -S ; }
  }

  // -------------------------------------------------------------------------
  // -------------------------------------------------------------------------

  void
  FresnelCS( indexType nk,
             valueType t,
             valueType C[],
             valueType S[] ) {
    FresnelCS(t,C[0],S[0]) ;
    if ( nk > 1 ) {
      valueType tt = m_pi_2*(t*t) ;
      valueType ss = sin(tt) ;
      valueType cc = cos(tt) ;
      C[1] = ss*m_1_pi ;
      S[1] = (1-cc)*m_1_pi ;
      if ( nk > 2 ) {
        C[2] = (t*ss-S[0])*m_1_pi ;
        S[2] = (C[0]-t*cc)*m_1_pi ;
      }
    }
  }

  //! \cond NODOC

  // -------------------------------------------------------------------------
  // -------------------------------------------------------------------------

  static
  void
  evalXYaLarge( valueType   a,
                valueType   b,
                valueType & X,
                valueType & Y ) {
    valueType s    = a > 0 ? +1 : -1 ;
    valueType absa = std::abs(a) ;
    valueType z    = m_1_sqrt_pi*sqrt(absa) ;
    valueType ell  = s*b*m_1_sqrt_pi/sqrt(absa) ;
    valueType g    = -0.5*s*(b*b)/absa ;
    valueType cg   = cos(g)/z ;
    valueType sg   = sin(g)/z ;

    valueType Cl, Sl, Cz, Sz ;
    FresnelCS( ell,   Cl, Sl ) ;
    FresnelCS( ell+z, Cz, Sz ) ;

    valueType dC0 = Cz - Cl ;
    valueType dS0 = Sz - Sl ;

    X = cg * dC0 - s * sg * dS0 ;
    Y = sg * dC0 + s * cg * dS0 ;
  }

  // -------------------------------------------------------------------------
  // nk max 3
  static
  void
  evalXYaLarge( indexType nk,
                valueType a,
                valueType b,
                valueType X[],
                valueType Y[] ) {

    G2LIB_ASSERT( nk < 4 && nk > 0,
                  "In evalXYaLarge first argument nk must be in 1..3, nk " << nk ) ;

    valueType s    = a > 0 ? +1 : -1 ;
    valueType absa = std::abs(a) ;
    valueType z    = m_1_sqrt_pi*sqrt(absa) ;
    valueType ell  = s*b*m_1_sqrt_pi/sqrt(absa) ;
    valueType g    = -0.5*s*(b*b)/absa ;
    valueType cg   = cos(g)/z ;
    valueType sg   = sin(g)/z ;

    valueType Cl[3], Sl[3], Cz[3], Sz[3] ;

    FresnelCS( nk, ell,   Cl, Sl ) ;
    FresnelCS( nk, ell+z, Cz, Sz ) ;

    valueType dC0 = Cz[0] - Cl[0] ;
    valueType dS0 = Sz[0] - Sl[0] ;
    X[0] = cg * dC0 - s * sg * dS0 ;
    Y[0] = sg * dC0 + s * cg * dS0 ;
    if ( nk > 1 ) {
      cg /= z ;
      sg /= z ;
      valueType dC1 = Cz[1] - Cl[1] ;
      valueType dS1 = Sz[1] - Sl[1] ;
      valueType DC  = dC1-ell*dC0 ;
      valueType DS  = dS1-ell*dS0 ;
      X[1] = cg * DC - s * sg * DS ;
      Y[1] = sg * DC + s * cg * DS ;
      if ( nk > 2 ) {
        valueType dC2 = Cz[2] - Cl[2] ;
        valueType dS2 = Sz[2] - Sl[2] ;
        DC   = dC2+ell*(ell*dC0-2*dC1) ;
        DS   = dS2+ell*(ell*dS0-2*dS1) ;
        cg   = cg/z ;
        sg   = sg/z ;
        X[2] = cg * DC - s * sg * DS ;
        Y[2] = sg * DC + s * cg * DS ;
      }
    }
  }

  // -------------------------------------------------------------------------
  // -------------------------------------------------------------------------

  valueType
  LommelReduced( valueType mu, valueType nu, valueType b ) {
    valueType tmp = 1/((mu+nu+1)*(mu-nu+1)) ;
    valueType res = tmp ;
    for ( indexType n = 1 ; n <= 100 ; ++n ) {
      tmp *= (-b/(2*n+mu-nu+1)) * (b/(2*n+mu+nu+1)) ;
      res += tmp ;
      if ( std::abs(tmp) < std::abs(res) * 1e-50 ) break ;
    }
    return res ;
  }

  // -------------------------------------------------------------------------
  // -------------------------------------------------------------------------

  static
  void
  evalXYazero( indexType nk,
               valueType b,
               valueType X[],
               valueType Y[] ) {

    valueType sb = sin(b) ;
    valueType cb = cos(b) ;
    valueType b2 = b*b ;
    if ( std::abs(b) < 1e-3 ) {
      X[0] = 1-(b2/6)*(1-(b2/20)*(1-(b2/42))) ;
      Y[0] = (b/2)*(1-(b2/12)*(1-(b2/30))) ;
    } else {
      X[0] = sb/b ;
      Y[0] = (1-cb)/b ;
    }
    // use recurrence in the stable part
    indexType m = indexType(floor(2*b)) ;
    if ( m >= nk ) m = nk-1 ;
    if ( m < 1   ) m = 1 ;
    for ( indexType k = 1 ; k < m ; ++k ) {
      X[k] = (sb-k*Y[k-1])/b ;
      Y[k] = (k*X[k-1]-cb)/b ;
    }
    //  use Lommel for the unstable part
    if ( m < nk ) {
      valueType A   = b*sb ;
      valueType D   = sb-b*cb ;
      valueType B   = b*D ;
      valueType C   = -b2*sb ;
      valueType rLa = LommelReduced(m+0.5,1.5,b) ;
      valueType rLd = LommelReduced(m+0.5,0.5,b) ;
      for ( indexType k = m ; k < nk ; ++k ) {
        valueType rLb = LommelReduced(k+1.5,0.5,b) ;
        valueType rLc = LommelReduced(k+1.5,1.5,b) ;
        X[k] = ( k*A*rLa + B*rLb + cb ) / (1+k) ;
        Y[k] = ( C*rLc + sb ) / (2+k) + D*rLd ;
	      rLa  = rLc ;
  	    rLd  = rLb ;
      }
    }
  }

  // -------------------------------------------------------------------------
  // -------------------------------------------------------------------------

  static
  void
  evalXYaSmall( valueType   a,
                valueType   b,
                indexType   p,
                valueType & X,
                valueType & Y ) {

    G2LIB_ASSERT( p < 11 && p > 0,
                  "In evalXYaSmall p = " << p << " must be in 1..10" ) ;

    valueType X0[43], Y0[43] ;

    indexType nkk = 4*p + 3 ; // max 43
    evalXYazero( nkk, b, X0, Y0 ) ;

    X = X0[0]-(a/2)*Y0[2] ;
    Y = Y0[0]+(a/2)*X0[2] ;

    valueType t  = 1 ;
    valueType aa = -a*a/4 ; // controllare!
    for ( indexType n=1 ; n <= p ; ++n ) {
      t *= aa/(2*n*(2*n-1)) ;
      valueType bf = a/(4*n+2) ;
      indexType jj = 4*n ;
      X += t*(X0[jj]-bf*Y0[jj+2]) ;
      Y += t*(Y0[jj]+bf*X0[jj+2]) ;
    }
  }

  // -------------------------------------------------------------------------

  static
  void
  evalXYaSmall( indexType nk,
                valueType a,
                valueType b,
                indexType p,
                valueType X[],
                valueType Y[] ) {

    indexType nkk = nk + 4*p + 2 ; // max 45
    valueType X0[45], Y0[45] ;

    G2LIB_ASSERT( nkk < 46,
                  "In evalXYaSmall (nk,p) = (" << nk << "," << p << ")\n" <<
                  "nk + 4*p + 2 = " << nkk  << " must be less than 46\n" ) ;

    evalXYazero( nkk, b, X0, Y0 ) ;

    for ( indexType j=0 ; j < nk ; ++j ) {
      X[j] = X0[j]-(a/2)*Y0[j+2] ;
      Y[j] = Y0[j]+(a/2)*X0[j+2] ;
    }

    valueType t  = 1 ;
    valueType aa = -a*a/4 ; // controllare!
    for ( indexType n=1 ; n <= p ; ++n ) {
      t *= aa/(2*n*(2*n-1)) ;
      valueType bf = a/(4*n+2) ;
      for ( indexType j = 0 ; j < nk ; ++j ) {
        indexType jj = 4*n+j ;
        X[j] += t*(X0[jj]-bf*Y0[jj+2]) ;
        Y[j] += t*(Y0[jj]+bf*X0[jj+2]) ;
      }
    }
  }
  
  //! \endcond

  // -------------------------------------------------------------------------
  // -------------------------------------------------------------------------

  void
  GeneralizedFresnelCS( valueType   a,
                        valueType   b,
                        valueType   c,
                        valueType & intC,
                        valueType & intS ) {

    valueType xx, yy ;
    if ( std::abs(a) < A_THRESOLD ) evalXYaSmall( a, b, A_SERIE_SIZE, xx, yy ) ;
    else                            evalXYaLarge( a, b, xx, yy ) ;

    valueType cosc = cos(c) ;
    valueType sinc = sin(c) ;

    intC = xx * cosc - yy * sinc ;
    intS = xx * sinc + yy * cosc ;
  }

  // -------------------------------------------------------------------------
  // -------------------------------------------------------------------------
  
  void
  GeneralizedFresnelCS( indexType nk,
                        valueType a,
                        valueType b,
                        valueType c,
                        valueType intC[],
                        valueType intS[] ) {

    G2LIB_ASSERT( nk > 0 && nk < 4, "nk = " << nk << " must be in 1..3" ) ;

    if ( std::abs(a) < A_THRESOLD ) evalXYaSmall( nk, a, b, A_SERIE_SIZE, intC, intS ) ;
    else                            evalXYaLarge( nk, a, b, intC, intS ) ;

    valueType cosc = cos(c) ;
    valueType sinc = sin(c) ;

    for ( indexType k = 0 ; k < nk ; ++k ) {
      valueType xx = intC[k] ;
      valueType yy = intS[k] ;
      intC[k] = xx * cosc - yy * sinc ;
      intS[k] = xx * sinc + yy * cosc ;
    }
  }

}

// EOF: Fresnel.cc
