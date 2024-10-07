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
/// file: Fresnel.hxx
///

namespace G2lib {

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /*\
   |   _____                         _
   |  |  ___| __ ___  ___ _ __   ___| |
   |  | |_ | '__/ _ \/ __| '_ \ / _ \ |
   |  |  _|| | |  __/\__ \ | | |  __/ |
   |  |_|  |_|  \___||___/_| |_|\___|_|
  \*/
  /*
   * Compute Fresnel integrals
   */
  void
  FresnelCS(
    real_type   x,
    real_type & C,
    real_type & S
  );

  //!
  //! Compute Fresnel integrals and its derivatives
  //!
  //! \f[
  //!   C(x) = \int_0^x \cos\left(\frac{\pi}{2}t^2\right) \mathrm{d}t, \qquad
  //!   S(x) = \int_0^x \sin\left(\frac{\pi}{2}t^2\right) \mathrm{d}t
  //! \f]
  //!
  //! \param nk maximum order of the derivative
  //! \param x  the input abscissa
  //! \param S  S[0]=\f$ S(x) \f$, S[1]=\f$ S'(x) \f$, S[2]=\f$ S''(x) \f$
  //! \param C  C[0]=\f$ C(x) \f$, C[1]=\f$ C'(x) \f$, C[2]=\f$ C''(x) \f$
  //!
  void
  FresnelCS(
    integer   nk,
    real_type x,
    real_type C[],
    real_type S[]
  );

  //!
  //! Compute the Fresnel integrals
  //!
  //! \f[
  //!   \int_0^1 t^k \cos\left(a\frac{t^2}{2} + b t + c\right)\,\mathrm{d}t,\qquad
  //!   \int_0^1 t^k \sin\left(a\frac{t^2}{2} + b t + c\right)\,\mathrm{d}t
  //! \f]
  //!
  //! \param nk   number of momentae to compute
  //! \param a    parameter \f$ a \f$
  //! \param b    parameter \f$ b \f$
  //! \param c    parameter \f$ c \f$
  //! \param intC cosine integrals,
  //! \param intS sine integrals
  //!
  void
  GeneralizedFresnelCS(
    integer   nk,
    real_type a,
    real_type b,
    real_type c,
    real_type intC[],
    real_type intS[]
  );

  //!
  //! Compute the Fresnel integrals
  //!
  //! \f[
  //!   \int_0^1 t^k \cos\left(a\frac{t^2}{2} + b t + c\right)\,\mathrm{d}t,\qquad
  //!   \int_0^1 t^k \sin\left(a\frac{t^2}{2} + b t + c\right)\,\mathrm{d}t
  //! \f]
  //!
  //! \param a      parameter \f$ a \f$
  //! \param b      parameter \f$ b \f$
  //! \param c      parameter \f$ c \f$
  //! \param intC   cosine integrals,
  //! \param intS   sine integrals
  //!
  void
  GeneralizedFresnelCS(
    real_type   a,
    real_type   b,
    real_type   c,
    real_type & intC,
    real_type & intS
  );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  //!
  //! Data storage for clothoid type curve
  //!
  class ClothoidData {
  public:

    real_type m_x0{0};     //!< initial \f$x\f$-coordinate of the clothoid
    real_type m_y0{0};     //!< initial \f$y\f$-coordinate of the clothoid
    real_type m_theta0{0}; //!< initial angle of the clothoid
    real_type m_kappa0{0}; //!< initial curvature
    real_type m_dk{0};     //!< curvature derivative

    ClothoidData() = default;

    real_type
    delta_theta( real_type s ) const
    { return s*(m_kappa0 + 0.5*s*m_dk); }

    real_type deltaTheta( real_type s ) const { return delta_theta( s ); }

    //!
    //! Return angle at curvilinear coordinate \f$s\f$
    //!
    real_type
    theta( real_type s ) const
    { return m_theta0 + s*(m_kappa0 + 0.5*s*m_dk); }

    real_type theta_D  ( real_type s ) const { return m_kappa0 + s*m_dk; }
    real_type theta_DD ( real_type   ) const { return m_dk; }
    real_type theta_DDD( real_type   ) const { return 0; }

    //!
    //! Return curvature at curvilinear coordinate \f$s\f$
    //!
    real_type kappa    ( real_type s ) const { return m_kappa0 + s*m_dk; }
    real_type kappa_D  ( real_type   ) const { return m_dk; }
    real_type kappa_DD ( real_type   ) const { return 0; }
    real_type kappa_DDD( real_type   ) const { return 0; }

    real_type X    ( real_type s ) const;
    real_type Y    ( real_type s ) const;
    real_type X_D  ( real_type s ) const;
    real_type Y_D  ( real_type s ) const;
    real_type X_DD ( real_type s ) const;
    real_type Y_DD ( real_type s ) const;
    real_type X_DDD( real_type s ) const;
    real_type Y_DDD( real_type s ) const;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    real_type X_ISO    ( real_type s, real_type offs ) const;
    real_type Y_ISO    ( real_type s, real_type offs ) const;
    real_type X_ISO_D  ( real_type s, real_type offs ) const;
    real_type Y_ISO_D  ( real_type s, real_type offs ) const;
    real_type X_ISO_DD ( real_type s, real_type offs ) const;
    real_type Y_ISO_DD ( real_type s, real_type offs ) const;
    real_type X_ISO_DDD( real_type s, real_type offs ) const;
    real_type Y_ISO_DDD( real_type s, real_type offs ) const;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    real_type X_SAE    ( real_type s, real_type offs ) const;
    real_type Y_SAE    ( real_type s, real_type offs ) const;
    real_type X_SAE_D  ( real_type s, real_type offs ) const;
    real_type Y_SAE_D  ( real_type s, real_type offs ) const;
    real_type X_SAE_DD ( real_type s, real_type offs ) const;
    real_type Y_SAE_DD ( real_type s, real_type offs ) const;
    real_type X_SAE_DDD( real_type s, real_type offs ) const;
    real_type Y_SAE_DDD( real_type s, real_type offs ) const;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    real_type tg0_x() const { return cos( m_theta0 ); }
    real_type tg0_y() const { return sin( m_theta0 ); }

    real_type tg_x( real_type s ) const { return cos(this->theta(s)); }
    real_type tg_y( real_type s ) const { return sin(this->theta(s)); }

    real_type tg_x_D( real_type s ) const;
    real_type tg_y_D( real_type s ) const;

    real_type tg_x_DD( real_type s ) const;
    real_type tg_y_DD( real_type s ) const;

    real_type tg_x_DDD( real_type s ) const;
    real_type tg_y_DDD( real_type s ) const;

    real_type nor0_x_ISO() const { return -this->tg0_y(); }
    real_type nor0_y_ISO() const { return this->tg0_x(); }

    real_type nor_x_ISO( real_type s ) const { return -this->tg_y(s); }
    real_type nor_y_ISO( real_type s ) const { return this->tg_x(s); }

    real_type nor_x_ISO_D( real_type s ) const { return -this->tg_y_D(s); }
    real_type nor_y_ISO_D( real_type s ) const { return this->tg_x_D(s);  }

    real_type nor_x_ISO_DD( real_type s ) const { return -this->tg_y_DD(s);  }
    real_type nor_y_ISO_DD( real_type s ) const { return this->tg_x_DD(s); }

    real_type nor_x_ISO_DDD( real_type s ) const { return -this->tg_y_DDD(s); }
    real_type nor_y_ISO_DDD( real_type s ) const { return this->tg_x_DDD(s); }

    real_type nor0_x_SAE() const { return this->tg0_y(); }
    real_type nor0_y_SAE() const { return -this->tg0_x(); }

    real_type nor_x_SAE( real_type s ) const { return this->tg_y(s); }
    real_type nor_y_SAE( real_type s ) const { return -this->tg_x(s); }

    real_type nor_x_SAE_D( real_type s ) const { return this->tg_y_D(s); }
    real_type nor_y_SAE_D( real_type s ) const { return -this->tg_x_D(s);  }

    real_type nor_x_SAE_DD( real_type s ) const { return this->tg_y_DD(s);  }
    real_type nor_y_SAE_DD( real_type s ) const { return -this->tg_x_DD(s); }

    real_type nor_x_SAE_DDD( real_type s ) const { return this->tg_y_DDD(s); }
    real_type nor_y_SAE_DDD( real_type s ) const { return -this->tg_x_DDD(s); }

    void tg( real_type s, real_type & tx, real_type & ty ) const;
    void tg_D( real_type s, real_type & tx, real_type & ty ) const;
    void tg_DD( real_type s, real_type & tx, real_type & ty ) const;
    void tg_DDD( real_type s, real_type & tx, real_type & ty ) const;

    void nor_ISO( real_type s, real_type & nx, real_type & ny ) const;
    void nor_ISO_D( real_type s, real_type & nx_D, real_type & ny_D ) const;
    void nor_ISO_DD( real_type s, real_type & nx_DD, real_type & ny_DD ) const;
    void nor_ISO_DDD( real_type s, real_type & nx_DDD, real_type & ny_DDD ) const;

    void nor_SAE( real_type s, real_type & nx, real_type & ny ) const;
    void nor_SAE_D( real_type s, real_type & nx_D, real_type & ny_D ) const;
    void nor_SAE_DD( real_type s, real_type & nx_DD, real_type & ny_DD ) const;
    void nor_SAE_DDD( real_type s, real_type & nx_DDD, real_type & ny_DDD ) const;

    void
    evaluate(
      real_type   s,
      real_type & theta,
      real_type & kappa,
      real_type & x,
      real_type & y
    ) const;

    void
    eval(
      real_type   s,
      real_type & x,
      real_type & y
    ) const;

    void
    eval_D(
      real_type   s,
      real_type & x_D,
      real_type & y_D
    ) const;

    void
    eval_DD(
      real_type   s,
      real_type & x_DD,
      real_type & y_DD
    ) const;

    void
    eval_DDD(
      real_type   s,
      real_type & x_DDD,
      real_type & y_DDD
    ) const;

    void
    eval_ISO(
      real_type   s,
      real_type   offs,
      real_type & x,
      real_type & y
    ) const;

    void
    eval_ISO_D(
      real_type   s,
      real_type   offs,
      real_type & x_D,
      real_type & y_D
    ) const;

    void
    eval_ISO_DD(
      real_type   s,
      real_type   offs,
      real_type & x_DD,
      real_type & y_DD
    ) const;

    void
    eval_ISO_DDD(
      real_type   s,
      real_type   offs,
      real_type & x_DDD,
      real_type & y_DDD
    ) const;

    void
    eval_SAE(
      real_type   s,
      real_type   offs,
      real_type & x,
      real_type & y
    ) const {
      this->eval_ISO( s, -offs, x, y );
    }

    void
    eval_SAE_D(
      real_type   s,
      real_type   offs,
      real_type & x_D,
      real_type & y_D
    ) const {
      this->eval_ISO_D( s, -offs, x_D, y_D );
    }

    void
    eval_DAE_DD(
      real_type   s,
      real_type   offs,
      real_type & x_DD,
      real_type & y_DD
    ) const {
      this->eval_ISO_DD( s, -offs, x_DD, y_DD );
    }

    void
    eval_SAE_DDD(
      real_type   s,
      real_type   offs,
      real_type & x_DDD,
      real_type & y_DDD
    ) const {
      this->eval_ISO_DDD( s, -offs, x_DDD, y_DDD );
    }

    void
    eval( real_type s, ClothoidData & C ) const;

    real_type c0x() const { return m_x0 - (sin(m_theta0)/m_kappa0); }
    real_type c0y() const { return m_y0 + (cos(m_theta0)/m_kappa0); }

    void
    Pinfinity( real_type & x, real_type & y, bool plus ) const;

    void
    reverse( real_type L );

    void
    reverse( real_type L, ClothoidData & out) const;

    void
    rotate( real_type angle, real_type cx, real_type cy );

    void
    origin_at( real_type s_origin );

    real_type
    split_at_flex( ClothoidData & C0, ClothoidData & C1 ) const;

    real_type
    aplus( real_type dtheta ) const;

    bool
    bbTriangle(
      real_type   L,
      real_type & xx0, real_type & yy0,
      real_type & xx1, real_type & yy1,
      real_type & xx2, real_type & yy2
    ) const;

    bool
    bbTriangle_ISO(
      real_type   L,
      real_type   offs,
      real_type & xx0, real_type & yy0,
      real_type & xx1, real_type & yy1,
      real_type & xx2, real_type & yy2
    ) const;

    bool
    bbTriangle_SAE(
      real_type   L,
      real_type   offs,
      real_type & xx0, real_type & yy0,
      real_type & xx1, real_type & yy1,
      real_type & xx2, real_type & yy2
    ) const {
      return this->bbTriangle_ISO( L, -offs, xx0, yy0, xx1, yy1, xx2, yy2 );
    }

    int
    build_G1(
      real_type   x0,
      real_type   y0,
      real_type   theta0,
      real_type   x1,
      real_type   y1,
      real_type   theta1,
      real_type   tol,
      real_type & L,
      bool        compute_deriv = false,
      real_type   L_D[2]        = nullptr,
      real_type   k_D[2]        = nullptr,
      real_type   dk_D[2]       = nullptr
    );

    bool
    build_forward(
      real_type   x0,
      real_type   y0,
      real_type   theta0,
      real_type   kappa0,
      real_type   x1,
      real_type   y1,
      real_type   tol,
      real_type & L
    );

    void
    info( ostream_type & s ) const;

  };

  #endif

}

///
/// eof: Fresnel.hxx
///
