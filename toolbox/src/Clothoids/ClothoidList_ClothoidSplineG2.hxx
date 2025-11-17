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
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

///
/// file: ClothoidList_ClothoidSplineG2.hxx
///

namespace G2lib {

  /*\
   |
   |    ___ _     _   _        _    _ ___      _ _           ___ ___
   |   / __| |___| |_| |_  ___(_)__| / __|_ __| (_)_ _  ___ / __|_  )
   |  | (__| / _ \  _| ' \/ _ \ / _` \__ \ '_ \ | | ' \/ -_) (_ |/ /
   |   \___|_\___/\__|_||_\___/_\__,_|___/ .__/_|_|_||_\___|\___/___|
   |                                     |_|
  \*/

  //!
  //! Class for the computation of \f$ G^2 \f$ spljne of clothoids
  //!
  class ClothoidSplineG2 {
  public:
    
    integer   m_max_iter{50};
    real_type m_tolerance{1e-10};
    real_type m_dump_min{0.1};

  private:

    TargetType m_tt{TargetType::P1};
    integer    m_npts{0};

    // work vector
    Eigen::Vector<real_type,Eigen::Dynamic> m_x;
    Eigen::Vector<real_type,Eigen::Dynamic> m_y;

    mutable vector<G2derivative> m_G2_vec;

    void evaluate_for_NLP( real_type const theta[] ) const;
    void evaluate_for_NLP_D( real_type const theta[] ) const;
    void evaluate_for_NLP_DD( real_type const theta[] ) const;

    void evaluate_for_NLP_BC( real_type const theta[] ) const;
    void evaluate_for_NLP_D_BC( real_type const theta[] ) const;
    void evaluate_for_NLP_DD_BC( real_type const theta[] ) const;

    real_type
    diff2pi( real_type in ) const
    { return in-Utils::m_2pi*round(in/Utils::m_2pi); }

    // vecchio da rimuovere
    void allocate( integer const npts );

  public:

    ClothoidSplineG2() = default;
    ~ClothoidSplineG2() = default;

    void setP4() { m_tt = TargetType::P4; }
    void setP5() { m_tt = TargetType::P5; }
    void setP6() { m_tt = TargetType::P6; }
    void setP7() { m_tt = TargetType::P7; }
    void setP8() { m_tt = TargetType::P8; }
    void setP9() { m_tt = TargetType::P9; }

    TargetType get_target() const { return m_tt; }

    // vecchio da rimuovere
    void
    build(
      real_type const xvec[],
      real_type const yvec[],
      integer   const npts
    );

    void
    build(
      integer   const npts,
      real_type const xvec[],
      real_type const yvec[],
      real_type       theta[]
    ) {
      build_PN( npts, xvec, yvec, theta, TargetType::P4 );
    }

    bool
    build_P1(
      integer   const npts,
      real_type const xvec[],
      real_type const yvec[],
      real_type       theta[],
      real_type       theta_init,
      real_type       theta_end
    );

    bool
    build_P2(
      integer   const npts,
      real_type const xvec[],
      real_type const yvec[],
      real_type       theta[]
    );

    bool
    build_PN(
      integer   const npts,
      real_type const xvec[],
      real_type const yvec[],
      real_type       theta[],
      TargetType      tt
    );

    integer numPnts() const { return m_npts; }
    integer numTheta() const;
    integer numConstraints() const;

    void
    guess(
      real_type theta_guess[],
      real_type theta_min[],
      real_type theta_max[]
    ) const;

    bool objective   ( real_type const theta[], real_type & f ) const;
    bool gradient    ( real_type const theta[], Pipal::Vector<real_type> & g ) const;
    bool constraints ( real_type const theta[], Pipal::Vector<real_type> & c ) const;

    bool jacobian ( real_type const theta[], Pipal::SparseMatrix<real_type> & J ) const;

    bool lagrangian_hessian ( real_type const theta[], real_type const lambda[], Pipal::SparseMatrix<real_type> & H ) const;

    string info() const;

    void info( ostream_type & stream ) const { stream << this->info(); }

    friend
    ostream_type &
    operator << ( ostream_type & stream, ClothoidSplineG2 const & c );

  };

}

///
/// eof: ClothoidList_ClothoidSplineG2.hxx
///
