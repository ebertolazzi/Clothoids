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

///
/// file: Dubins3p_ellipse.cc
///

#include "Clothoids.hh"
#include "Clothoids_fmt.hh"

#include "PolynomialRoots.hh"

namespace G2lib {

  using PolynomialRoots::Quartic;

  bool
  Dubins3p::build_ellipse(
    real_type xi,
    real_type yi,
    real_type thetai,
    real_type xm,
    real_type ym,
    real_type xf,
    real_type yf,
    real_type thetaf,
    real_type k_max
  ) {

    m_evaluation = 0;

    typedef struct Dubins3p_data {
      Dubins    D0{"temporary Dubins A"};
      Dubins    D1{"temporary Dubins B"};
      real_type thetam{0};
      real_type len{0};

      void
      copy( Dubins3p_data const & rhs ) {
        D0.copy(rhs.D0);
        D1.copy(rhs.D1);
        thetam = rhs.thetam;
        len    = rhs.len;
      }

      //bool
      //compare( Dubins3p_data const & D ) const {
      //  return D0.solution_type() == D.D0.solution_type() &&
      //         D1.solution_type() == D.D1.solution_type();
      //}

    } Dubins3p_data;

    auto eval3p = [this,xi,yi,thetai,xm,ym,xf,yf,thetaf,k_max]( Dubins3p_data & D3P ) -> void {
      D3P.D0.build( xi, yi, thetai, xm, ym, D3P.thetam, k_max );
      D3P.D1.build( xm, ym, D3P.thetam, xf, yf, thetaf, k_max );
      D3P.len = D3P.D0.length() + D3P.D1.length();
      ++m_evaluation;
    };

    real_type min_len{ Utils::Inf<real_type>() };

    real_type const RR{1/k_max};

    real_type const ii[4]{RR, RR,-RR,-RR};
    real_type const jj[4]{RR,-RR, RR,-RR};

    // rotate -90 degree the tangent
    real_type const ni[2]{sin(thetai),-cos(thetai)};
    real_type const nf[2]{sin(thetaf),-cos(thetaf)};

    Dubins3p_data D3P;

    for ( integer kk{0}; kk < 4; ++kk ) {

      // centro cerchi
      real_type const dI{ii[kk]};
      real_type const dF{jj[kk]};
      real_type const cxI{xi+dI*ni[0]};
      real_type const cyI{yi+dI*ni[1]};
      real_type const cxF{xf+dF*nf[0]};
      real_type const cyF{yf+dF*nf[1]};

      // setup in standard form (-1,0), (1,0)
      real_type const dx{ cxF - cxI };
      real_type const dy{ cyF - cyI };
      real_type const th{ atan2( dy, dx ) };
      real_type const l { hypot( dx, dy ) };
      real_type const l2{ l*l };

      real_type const X{ xm - cxI };
      real_type const Y{ ym - cyI };

      // proietta e ruota
      real_type const x0 { (2/l2)*(X*dx+Y*dy)-1 };
      real_type const y0 { (2/l2)*(Y*dx-X*dy)   };
      real_type const R  { (2/l)*RR };

      // calcolo polinomio per cercare candidati
      real_type const x02{x0*x0};
      real_type const y02{y0*y0};
      real_type const A{x0*(R-y0)};
      real_type const B{-2*R*y0 - 2*x02 + 2*y02 + 2};
      real_type const C{6*x0*y0};
      real_type const D{-2*R*y0 + 2*x02 - 2*y02 - 2};
      real_type const E{-x0*(R+y0)};

      Quartic Q(A,B,C,D,E);

      real_type     r[4];
      integer const nr{ Q.get_real_roots( r ) };
      for ( integer i{0}; i < nr; ++i ) {

        real_type const theta(2*atan(r[i]));
        real_type const cc{ cos(theta) };
        real_type const ss{ sin(theta) };
        real_type a{ -(R*cc+y0)*ss/(ss*y0+cc*x0) };

        if ( a < 0 ) continue;

        // ritorna alla soluzione principale
        D3P.thetam = theta + th;

        // calcolo Dubins composita
        eval3p( D3P );

        if ( D3P.len < min_len ) {
          min_len = D3P.len;
          m_Dubins0.copy(D3P.D0);
          m_Dubins1.copy(D3P.D1);
        }
      }
    }

    return min_len != Utils::Inf<real_type>();
  }

}

///
/// eof: Dubins3p_ellipse.cc
///
