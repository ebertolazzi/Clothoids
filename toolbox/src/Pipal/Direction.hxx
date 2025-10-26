/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Pipal project is distributed under the MIT License.                                       *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef INCLUDE_PIPAL_DIRECTION_HXX
#define INCLUDE_PIPAL_DIRECTION_HXX

namespace Pipal {

  /**
   * \brief Reset a search direction to zero and initialize norms.
   * \tparam Real Floating-point type used by the algorithm.
   * \param[in] d Direction object to reset.
   * \param[in] i Input structure providing vector/matrix sizes.
   */
  template<typename Real>
  inline
  void
  Solver<Real>::resetDirection( Direction<Real> & d ) const {
    Input<Real> const & i{this->m_input};
    d.x.setZero(i.nV);
    d.r1.setZero(i.nE);
    d.r2.setZero(i.nE);
    d.lE.setZero(i.nE);
    d.s1.setZero(i.nI);
    d.s2.setZero(i.nI);
    d.lI.setZero(i.nI);
    d.x_norm  = 0;
    d.x_norm_ = std::numeric_limits<Real>::infinity();
    d.l_norm  = 0;
    d.lred0 = d.ltred0 = d.ltred = d.qtred = d.m = 0;
  }

  /**
   * \brief Evaluate a linear combination of up to three directions.
   * \tparam Real Floating-point type used by the algorithm.
   * \param[in] d Output direction where the result is stored.
   * \param[in] i Input structure (used for size checks).
   * \param[in] d1 First direction.
   * \param[in] d2 Second direction.
   * \param[in] d3 Third direction.
   * \param[in] a1 Scalar multiplier for the first direction.
   * \param[in] a2 Scalar multiplier for the second direction.
   * \param[in] a3 Scalar multiplier for the third direction.
   */
  template<typename Real>
  inline
  void
  Solver<Real>::evalLinearCombination(
    Direction<Real> const & d1,
    Direction<Real> const & d2,
    Direction<Real> const & d3,
    Real            const   a1,
    Real            const   a2,
    Real            const   a3
  ) {
    // Create alias for easier access
    Input<Real>     & i{this->m_input};
    Direction<Real> & d{this->m_direction};

    // Evaluate linear combinations
    d.x = a1*d1.x + a2*d2.x + a3*d3.x;
    if (i.nE > 0) {
      d.r1 = a1*d1.r1 + a2*d2.r1 + a3*d3.r1;
      d.r2 = a1*d1.r2 + a2*d2.r2 + a3*d3.r2;
      d.lE = a1*d1.lE + a2*d2.lE + a3*d3.lE;
    }
    if (i.nI > 0) {
      d.s1 = a1*d1.s1 + a2*d2.s1 + a3*d3.s1;
      d.s2 = a1*d1.s2 + a2*d2.s2 + a3*d3.s2;
      d.lI = a1*d1.lI + a2*d2.lI + a3*d3.lI;
    }

    // Evaluate primal direction norm
    d.x_norm = d.x.norm();

    // Evaluate dual direction norm
    d.l_norm = std::sqrt(d.lE.matrix().squaredNorm() + d.lI.matrix().squaredNorm());
  }

  /**
   * \brief Evaluate model reductions and quality metric for a direction.
   * \tparam Real Floating-point type used by the algorithm.
   * \param[in] d Direction to analyze (most fields updated).
   * \param[in] i Input structure containing sizes.
   * \param[in] z Current iterate with model data (gradients, Jacobians, etc.).
   */
  template<typename Real>
  inline
  void
  Solver<Real>::evalModels() {
    // Create alias for easier access
    Input<Real>     & i{this->m_input};
    Iterate<Real>   & z{this->m_iterate};
    Direction<Real> & d{this->m_direction};
  
    // Evaluate reduction in linear model of penalty-interior-point objective for zero penalty parameter
    d.lred0 = 0;
    if (i.nE > 0) {
      d.lred0 -= ((1.0 - z.mu / z.r1) * d.r1).sum();
      d.lred0 -= ((1.0 - z.mu / z.r2) * d.r2).sum();
    }
    if (i.nI > 0) {
      d.lred0 -= (((-z.mu) / z.s1) * d.s1).sum();
      d.lred0 -= ((1.0 - z.mu / z.s2) * d.s2).sum();
    }

    // Evaluate remaining quantities only for nonzero penalty parameter
    if (z.rho > 0)
    {
      // Evaluate reduction in linear model of merit function for zero penalty parameter
      d.ltred0 = 0;
      if (i.nE > 0) {
        Array<Real> sqrt_term((z.cE.square() + z.mu*z.mu).sqrt());
        d.ltred0 -= 0.5*(((1.0 - z.mu/z.r1) * (-1.0 + z.cE/sqrt_term) + (1.0 - z.mu/z.r2) *
          (1.0 + z.cE/sqrt_term)) * (z.JE * d.x).array()).sum();
      }
      if (i.nI > 0) {
        Array<Real> sqrt_term((z.cI.square() + 4.0*z.mu*z.mu).sqrt());
        d.ltred0 -= 0.5*((((-z.mu)/z.s1) * (-1.0 + z.cI/sqrt_term) + (1.0 - z.mu/z.s2) *
          (1.0 + z.cI/sqrt_term)) * (z.JI * d.x).array()).sum();
      }

      // Evaluate reduction in linear model of merit function
      d.ltred = -z.rho*z.g.transpose()*d.x + d.ltred0;

      // Evaluate reduction in quadratic model of merit function
      d.qtred = d.ltred - 0.5*d.x.transpose()*z.H*d.x;
      if (i.nE > 0) {
        Array<Real> Jd(z.JE*d.x);
        Array<Real> Dinv((z.r1/(1.0+z.lE) + z.r2/(1.0-z.lE)).matrix());
        d.qtred -= 0.5*Jd.matrix().transpose() * ((Jd/Dinv).matrix());
      }
      if (i.nI > 0) {
        Array<Real> Jd(z.JI*d.x);
        Array<Real> Dinv((z.s1/(0.0+z.lI) + z.s2/(1.0-z.lI)).matrix());
        d.qtred -= 0.5*Jd.matrix().transpose() * ((Jd/Dinv).matrix());
      }

      // Initialize quality function vector
      Array<Real> vec(i.nV+2*i.nE+2*i.nI);
      vec.setZero();

      // Set gradient of objective
      vec.head(i.nV) = z.rho*z.g;

      // Set gradient of Lagrangian for constraints
      if (i.nE > 0) {vec.head(i.nV) += ((z.lE+d.lE).matrix().transpose()*z.JE).transpose().array();}
      if (i.nI > 0) {vec.head(i.nV) += ((z.lI+d.lI).matrix().transpose()*z.JI).transpose().array();}

      // Set complementarity for constraint slacks
      if (i.nE > 0) {
          vec.segment(i.nV, 2*i.nE) <<
            (z.r1+d.r1)*(1.0 + (z.lE+d.lE)), (z.r2+d.r2)*(1.0 - (z.lE+d.lE));
      }
      if (i.nI > 0) {
          vec.segment(i.nV+2*i.nE, 2*i.nI) <<
            (z.s1+d.s1)*(z.lI+d.lI), (z.s2+d.s2)*(1.0 - (z.lI+d.lI));
      }

      // Evaluate quality function
      d.m = vec.matrix().template lpNorm<Eigen::Infinity>();
    }
  }

  /**
   * \brief Recover direction components from the Newton system solution.
   * \tparam Real Floating-point type used by the algorithm.
   * \param[in] d Direction object whose components are set from the Newton solution.
   */
  template<typename Real>
  inline
  void
  Solver<Real>::evalNewtonStep() {
    // Create alias for easier access
    Input<Real>   const & i{this->m_input};
    Iterate<Real> const & z{this->m_iterate};
    Direction<Real>     & d{this->m_direction};

    // Evaluate direction
    SparseVector<Real> dir(z.ldlt.solve(-z.b));

    // Parse direction
    d.x = dir.head(i.nV);
    if (i.nE > 0) {
      d.r1 = dir.segment(i.nV, i.nE);
      d.r2 = dir.segment(i.nV+i.nE, i.nE);
      d.lE = dir.segment(i.nV+i.nE+i.nE+i.nI+i.nI, i.nE);
    }
    if (i.nI > 0) {
      d.s1 = dir.segment(i.nV+i.nE+i.nE, i.nI);
      d.s2 = dir.segment(i.nV+i.nE+i.nE+i.nI, i.nI);
      d.lI = dir.segment(i.nV+i.nE+i.nE+i.nI+i.nI+i.nE, i.nI);
    }

    // Evaluate primal direction norm
    d.x_norm = d.x.norm();

    // Evaluate dual direction norm
    d.l_norm = std::sqrt(d.lE.matrix().squaredNorm() + d.lI.matrix().squaredNorm());
  }

  /**
   * \brief Compute the search direction for the current iterate.
   * \tparam Real Floating-point type used by the algorithm.
   * \param[in] d Direction to fill.
   * \param[in] p Algorithm parameters (may be updated).
   * \param[in] i Problem input structure.
   * \param[in] c Counters used for bookkeeping of evaluations.
   * \param[in] z Current iterate and working matrices.
   * \param[in] a Acceptance object used by some subroutines.
   * \param[in] problem Problem interface used for evaluations.
   */
  template <typename Real>
  inline
  void
  Solver<Real>::evalStep() {

    // Create alias for easier access
    Parameter<Real> & p{this->m_parameter};
    Iterate<Real>   & z{this->m_iterate};
    Direction<Real> & d{this->m_direction};

    // Reset maximum exponent for interior-point parameter increases
    resetMuMaxExp();

    // Update penalty-interior-point parameters based on KKT errors
    updateParameters();

    // Evaluate matrices
    evalMatrices();

    // Set last penalty parameter
    setRhoLast(z.rho);

    // Check for aggressive algorithm
    if (p.algorithm == Algorithm::ADAPTIVE)
    {
      // Check KKT memory for potential mu increase limit
      if (z.kkt(1) > z.kkt_.maxCoeff()) setMuMaxExpZero();

      // Store current penalty and interior-point parameters
      Real rho_curr{z.rho}, mu_curr{z.mu};

      // Evaluate trial steps
      Direction<Real> d1, d2, d3;
      resetDirection(d1);
      resetDirection(d2);
      resetDirection(d3);
      evalTrialSteps(d1, d2, d3);

      // Set trial interior-point parameter values
      Array<Real> exponents(Array<Real>::LinSpaced(p.mu_trials, p.mu_trials - 1, 0) - p.mu_max_exp);
      Array<Real> Mu((mu_curr * exponents.unaryExpr(
        [&p] (Real e) {return std::pow(p.mu_factor, e);})
        ).min(p.mu_max).max(p.mu_min));

      // Initialize feasibility direction data
      Vector<Real> lred0_0_mu(p.mu_trials);
      lred0_0_mu.setZero();

      // Loop through interior-point parameter values
      for (Integer j{0}; j < p.mu_trials; ++j)
      {
        // Set penalty and interior-point parameters
        setRho(0.0);
        setMu(Mu(j));

        // Evaluate direction
        evalLinearCombination(
          d1, d2, d3,
          z.rho/rho_curr+z.mu/mu_curr-1.0,
          1.0-z.mu/mu_curr,
          1.0-z.rho/rho_curr
        );

        // Cut length
        d.x *= std::min(d.x_norm_/std::max(d.x_norm, 1.0), 1.0);

        // Run fraction-to-boundary
        fractionToBoundary();

        // Cut length
        evalTrialStepCut();

        // Evaluate models
        evalModels();

        // Set feasibility direction data
        lred0_0_mu(j) = d.lred0;
      }

      // Initialize updating data
      Vector<Real> ltred0_rho_mu(p.mu_trials), qtred_rho_mu(p.mu_trials), m_rho_mu(p.mu_trials);
      ltred0_rho_mu.setZero(); qtred_rho_mu.setZero(); m_rho_mu.setZero();

      // Initialize check
      bool check{false};

      // Loop through penalty parameter values
      for (Integer k{0}; k < p.rho_trials; ++k)
      {
        // Set penalty parameter
        setRho( std::max(p.rho_min, std::pow(p.rho_factor, k)*rho_curr) );

        // Set last penalty parameter
        if (rho_curr > z.kkt(0)*z.kkt(0)) setRhoLast( z.rho );

        // Loop through interior-point parameter values
        for (Integer j{0}; j < p.mu_trials; ++j)
        {
          // Set interior-point parameter
          setMu(Mu(j));

          // Evaluate direction
          evalLinearCombination(
            d1, d2, d3,
            z.rho/rho_curr+z.mu/mu_curr-1.0,
            1.0-z.mu/mu_curr,
            1.0-z.rho/rho_curr
          );

          // Run fraction-to-boundary
          fractionToBoundary();

          // Cut steps
          evalTrialStepCut();

          // Evaluate models
          evalModels();

          // Set updating data
          ltred0_rho_mu(j) = d.ltred0;
          qtred_rho_mu(j)  = d.qtred;
          m_rho_mu(j)      = d.m;

          // Check updating conditions for infeasible points
          if (z.v > p.opt_err_tol && (ltred0_rho_mu(j) < p.update_con_1*lred0_0_mu(j) ||
            qtred_rho_mu(j) < p.update_con_2*lred0_0_mu(j) || z.rho > z.kkt(0)*z.kkt(0))) {
            m_rho_mu(j) = std::numeric_limits<Real>::infinity();
          }

          // Check updating conditions for feasible points
          if (z.v <= p.opt_err_tol && qtred_rho_mu(j) < 0.0) {
            m_rho_mu(j) = std::numeric_limits<Real>::infinity();
          }
        }

        // Find minimum m for current rho
        Real m_min{m_rho_mu.minCoeff()};

        // Check for finite minimum
        if (m_min < std::numeric_limits<Real>::infinity())
        {
          // Loop through mu values
          for (Integer j{0}; j < p.mu_trials; ++j)
          {
            // Check condition
            if (m_rho_mu(j) <= p.update_con_3*m_min) setMu(Mu(j));
          }

          // Set condition check
          check = true;

          // Break loop
          break;
        }
      }

      // Check conditions
      if (check == false) { setRho(rho_curr); setMu(mu_curr); }

      // Evaluate merit
      evalMerit();
    }

    // Evaluate primal-dual right-hand side vector
    evalNewtonRhs();

    // Evaluate search direction
    evalNewtonStep();

    // Evaluate models
    evalModels();

    // Store last direction norm
    d.x_norm_ = d.x_norm;
  }

  /**
   * \brief Store a trial step into another direction object.
   * \tparam Real Floating-point type used by the algorithm.
   * \param[out] v Destination direction that will contain the trial step.
   */
  template<typename Real>
  inline
  void
  Solver<Real>::evalTrialStep( Direction<Real> & v ) const {
    Input<Real>     const & i{this->m_input};
    Direction<Real> const & d{this->m_direction};
    // Set direction components
    v.x = d.x;
    if (i.nE > 0) {v.r1 = d.r1; v.r2 = d.r2; v.lE = d.lE;}
    if (i.nI > 0) {v.s1 = d.s1; v.s2 = d.s2; v.lI = d.lI;}
  }

  /**
   * \brief Scale a trial step by the fraction-to-boundary values.
   * \tparam Real Floating-point type used by the algorithm.
   * \param[in] d Direction to scale.
   * \param[in] i Input structure (used for sizes).
   * \param[in] a Acceptance object containing fraction values (p and d).
   */
  template<typename Real>
  inline
  void
  Solver<Real>::evalTrialStepCut() {
    // Create alias for easier access
    Input<Real>      & i{this->m_input};
    Direction<Real>  & d{this->m_direction};
    Acceptance<Real> & a{this->m_acceptance};

    // Set direction components
    d.x = a.p*d.x ;
    if (i.nE > 0) {d.r1 = a.p*d.r1; d.r2 = a.p*d.r2; d.lE = a.d*d.lE;}
    if (i.nI > 0) {d.s1 = a.p*d.s1; d.s2 = a.p*d.s2; d.lI = a.d*d.lI;}
  }

  /**
   * \brief Compute and store directions for a few parameter combinations.
   * \tparam Real Floating-point type used by the algorithm.
   * \param[out] d1 Direction for the first parameter combination.
   * \param[out] d2 Direction for the second parameter combination.
   * \param[out] d3 Direction for the third parameter combination.
   */
  template<typename Real>
  inline
  void
  Solver<Real>::evalTrialSteps(
    Direction<Real> & d1,
    Direction<Real> & d2,
    Direction<Real> & d3
  ) {
    // Create alias for easier access
    Iterate<Real> const & z{this->m_iterate};

    // Store current penalty and interior-point parameters
    Real rho_curr{z.rho}, mu_curr{z.mu};

    // Evaluate direction for current penalty and interior-point parameters
    setRho(rho_curr);
    setMu(mu_curr);
    evalNewtonRhs();
    evalNewtonStep();
    evalTrialStep(d1);

    // Evaluate direction for zero interior-point parameter
    setRho(rho_curr);
    setMu(0.0);
    evalNewtonRhs();
    evalNewtonStep();
    evalTrialStep(d2);

    // Evaluate direction for zero penalty parameter
    setRho(0.0);
    setMu(mu_curr);
    evalNewtonRhs();
    evalNewtonStep();
    evalTrialStep(d3);
  }

  /**
   * \brief Populate a Direction object from its component vectors.
   * \tparam Real Floating-point type used by the algorithm.
   * \param[out] d Direction to populate.
   * \param[in] i Input used to check which blocks are present.
   * \param[in] dx Primal component of the direction.
   * \param[in] dr1 Equality slack direction (first part).
   * \param[in] dr2 Equality slack direction (second part).
   * \param[in] ds1 Inequality slack direction (first part).
   * \param[in] ds2 Inequality slack direction (second part).
   * \param[in] dlE Equality multipliers direction.
   * \param[in] dlI Inequality multipliers direction.
   * \param[in] dx_norm Norm of the primal direction (precomputed).
   * \param[in] dl_norm Norm of the dual direction (precomputed).
   */
  template<typename Real>
  inline
  void
  Solver<Real>::setDirection(
    Vector<Real> const & dx,
    Vector<Real> const & dr1,
    Vector<Real> const & dr2,
    Vector<Real> const & ds1,
    Vector<Real> const & ds2,
    Vector<Real> const & dlE,
    Vector<Real> const & dlI,
    Real         const   dx_norm,
    Real         const   dl_norm
  ) {
    // Create alias for easier access
    Input<Real> const & i{this->m_input};
    Direction<Real>   & d{this->m_direction};

    // Set primal variables
    d.x = dx;
    if (i.nE > 0) {d.r1 = dr1; d.r2 = dr2; d.lE = dlE;}
    if (i.nI > 0) {d.s1 = ds1; d.s2 = ds2; d.lI = dlI;}
    d.x_norm = dx_norm;
    d.l_norm = dl_norm;
  }

} // namespace Pipal

#endif /* INCLUDE_PIPAL_DIRECTION_HH */
