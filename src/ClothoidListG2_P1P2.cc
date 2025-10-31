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

#include "Clothoids.hh"
#include "Clothoids_fmt.hh"

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wswitch-enum"
#endif
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wswitch-enum"
#endif

namespace G2lib {

  using Vector       = Pipal::Vector<real_type>;
  using SparseMatrix = Pipal::SparseMatrix<real_type>;

  /*\
   |
   |    ___ _     _   _        _    _ ___      _ _           ___ ___
   |   / __| |___| |_| |_  ___(_)__| / __|_ __| (_)_ _  ___ / __|_  )
   |  | (__| / _ \  _| ' \/ _ \ / _` \__ \ '_ \ | | ' \/ -_) (_ |/ /
   |   \___|_\___/\__|_||_\___/_\__,_|___/ .__/_|_|_||_\___|\___/___|
   |                                     |_|
  \*/
  
  static auto Newton = [](
    integer   n,
    integer   max_iter,
    real_type dump_min,
    real_type tolerance,
    real_type theta[],
    auto &    F,
    auto &    F_JF
  ) -> bool {
    // applico Newton dumped
    Vector x0( n ), x1( n ), f0( n ), f1( n ), d0( n ), d1( n );
    std::copy_n( theta, n, x0.data() );
    bool converged{false};
    Eigen::SparseLU<SparseMatrix> solver;
    for ( integer iter{0}; iter < max_iter && !converged; ++iter ) {
      SparseMatrix JF0;
      if ( !F_JF( x0, f0, JF0 ) ) break;
      // calcola direzione di ricerca
      solver.analyzePattern(JF0);  // Fase simbolica (struttura)
      solver.factorize(JF0);       // Fase numerica (valori)
      d0 = solver.solve(f0);
      real_type norm_d0{ d0.norm() / n };
      real_type norm_f0{ f0.lpNorm<1>() / n };
      real_type dump_reduction{0.75};
      real_type alpha{1/dump_reduction};
      do {
        alpha *= dump_reduction;
        x1.noalias() = x0 - alpha * d0;
        // devo aggiustare gli angoli in modo che differiscano meno di 2*pi
        //for ( integer i{1}; i < m_npts; ++i ) adjust( x1(i-1), x1(i) );
        bool ok{ F( x1, f1 ) };
        if ( ok ) {
          d1 = solver.solve(f1);
          real_type norm_d1{ d1.norm() / n };
          real_type norm_f1{ f1.lpNorm<1>() / n };
          if ( norm_d1 <= norm_d0*(1-alpha/2) || norm_f1 < 0.99 * norm_f0 ) break;
        }
      } while ( alpha > dump_min );
      // fmt::print( "it:{:2} α={} ‖f‖:{}\n", iter, alpha, norm_f0 );
      x0 = x1;
      converged = norm_f0 < tolerance;
    }
    std::copy_n( x0.data(), n, theta );
    return converged;
  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidSplineG2::build_P1(
    integer   const n,
    real_type const xvec[],
    real_type const yvec[],
    real_type       theta[],
    real_type       theta_init,
    real_type       theta_end
  ) {
  
    allocate( n );

    std::copy_n( xvec, m_npts, m_x.data() );
    std::copy_n( yvec, m_npts, m_y.data() );

    {
      Vector theta_min( m_npts ),
             theta_max( m_npts ),
             omega( m_npts ),
             len( m_npts );
      G2lib::xy_to_guess_angle( m_npts, xvec, yvec, theta, theta_min.data(), theta_max.data(), omega.data(), len.data() );
    }

    integer const ne  { m_npts - 1 };
    integer const ne1 { m_npts - 2 };

    auto F = [this,&ne,&ne1,&theta_init,&theta_end]( Vector const & theta, Vector & F ) -> bool {
      this->evaluate_for_NLP( theta.data() );
      for ( integer j{0}; j < ne1; ++j ) F[j] = m_G1_vec[j].k1-m_G1_vec[j+1].k0;
      F[ne1] = diff2pi( theta[0]  - theta_init );
      F[ne]  = diff2pi( theta[ne] - theta_end );
      return F.allFinite();;
    };

    auto F_JF = [this,&ne,&ne1,&theta_init,&theta_end]( Vector const & theta, Vector & F, SparseMatrix & JF ) -> bool {
      JF.setZero();
      JF.resize( m_npts, m_npts );

      this->evaluate_for_NLP_D( theta.data() );
      for ( integer j{0}; j < ne1; ++j ) F[j] = m_G1_vec[j].k1 - m_G1_vec[j+1].k0;
      F[ne1] = diff2pi( theta[0]  - theta_init );
      F[ne]  = diff2pi( theta[ne] - theta_end );

      vector<Eigen::Triplet<real_type>> triplets;
      triplets.reserve(3*ne+3);
      for ( integer j{0}; j < ne1; ++j ) {
        auto const & G  { m_G1_vec[j] };
        auto const & G1 { m_G1_vec[j+1] };
        triplets.emplace_back(j, j,   G.k__L + G.dk__L * G.L + G.dk * G.L__L           );
        triplets.emplace_back(j, j+1, G.k__R + G.dk__R * G.L + G.dk * G.L__R - G1.k__L );
        triplets.emplace_back(j, j+2,                                        - G1.k__R );
      }
      triplets.emplace_back( ne1, 0,  1 );
      triplets.emplace_back( ne,  ne, 1 );

      JF.setFromTriplets(triplets.begin(), triplets.end());
      JF.makeCompressed();

      return F.allFinite();
    };
    
    return Newton( n, m_max_iter, m_dump_min, m_tolerance, theta, F, F_JF );

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidSplineG2::build_P2(
    integer   const n,
    real_type const xvec[],
    real_type const yvec[],
    real_type       theta[]
  ) {
  
    allocate( n );

    std::copy_n( xvec, m_npts, m_x.data() );
    std::copy_n( yvec, m_npts, m_y.data() );

    {
      Vector theta_min( m_npts ),
             theta_max( m_npts ),
             omega( m_npts ),
             len( m_npts );
      G2lib::xy_to_guess_angle( m_npts, xvec, yvec, theta, theta_min.data(), theta_max.data(), omega.data(), len.data() );
    }

    integer const ne  { m_npts - 1 };
    integer const ne1 { m_npts - 2 };

    auto F = [this,&ne,&ne1]( Vector const & theta, Vector & F ) -> bool {
      this->evaluate_for_NLP( theta.data() );
      for ( integer j{0}; j < ne1; ++j ) F[j] = m_G1_vec[j].k1-m_G1_vec[j+1].k0;
      F[ne1] = diff2pi( theta[0] - theta[ne] );
      F[ne]  = m_G1_vec[ne1].k1-m_G1_vec[0].k0;
      return F.allFinite();;
    };

    auto F_JF = [this,&ne,&ne1]( Vector const & theta, Vector & F, SparseMatrix & JF ) -> bool {
      JF.setZero();
      JF.resize( m_npts, m_npts );

      this->evaluate_for_NLP_D( theta.data() );
      for ( integer j{0}; j < ne1; ++j ) F[j] = m_G1_vec[j].k1 - m_G1_vec[j+1].k0;
      F[ne1] = diff2pi( theta[0] - theta[ne] );
      F[ne]  = m_G1_vec[ne1].k1-m_G1_vec[0].k0;

      vector<Eigen::Triplet<real_type>> triplets;
      triplets.reserve(3*ne+4);
      for ( integer j{0}; j < ne1; ++j ) {
        auto const & G  { m_G1_vec[j] };
        auto const & G1 { m_G1_vec[j+1] };
        triplets.emplace_back(j, j,   G.k__L + G.dk__L * G.L + G.dk * G.L__L           );
        triplets.emplace_back(j, j+1, G.k__R + G.dk__R * G.L + G.dk * G.L__R - G1.k__L );
        triplets.emplace_back(j, j+2,                                        - G1.k__R );
      }
      triplets.emplace_back( ne1, 0,   1 );
      triplets.emplace_back( ne1, ne, -1 );
      {
        auto const & G  { m_G1_vec[ne1] };
        auto const & G1 { m_G1_vec[0]   };
        triplets.emplace_back(ne, ne1, G.k__L + G.dk__L * G.L + G.dk * G.L__L );
        triplets.emplace_back(ne, ne,  G.k__R + G.dk__R * G.L + G.dk * G.L__R );
        triplets.emplace_back(ne, 0,   - G1.k__L );
        triplets.emplace_back(ne, 1,   - G1.k__R );
      }

      JF.setFromTriplets(triplets.begin(), triplets.end());
      JF.makeCompressed();

      return F.allFinite();
    };

    return Newton( n, m_max_iter, m_dump_min, m_tolerance, theta, F, F_JF );

  }

}

// EOF: ClothoidListG2_P1P2.cc
