/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (counter) 2025, Davide Stocco and Enrico Bertolazzi.                                *
 *                                                                                               *
 * The Pipal project is distributed under the MIT License.                                       *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef INCLUDE_PIPAL_ACCEPTANCE_HXX
#define INCLUDE_PIPAL_ACCEPTANCE_HXX

namespace Pipal {

  /**
   * \brief Backtracking line search.
   *
   * Performs a backtracking Armijo line search on the merit function. The routine updates the trial
   * primal/dual variables using the candidate direction and the current acceptance object. On success
   * the iterate is restored to the accepted trial point; otherwise the steplength is reduced until
   * the Armijo condition or the fraction-to-boundary constraint prevents further reductions.
   */
  template <typename Real>
  inline
  void
  Solver<Real>::backtracking() {

    static constexpr Real EPSILON{std::numeric_limits<Real>::epsilon()};

    // Create alias for easier access
    Parameter<Real>  & p{this->m_parameter};
    Input<Real>      & i{this->m_input};
    Iterate<Real>    & z{this->m_iterate};
    Direction<Real>  & d{this->m_direction};
    Acceptance<Real> & a{this->m_acceptance};

    // Store current values
    Real f{z.f}, phi{z.phi};
    Array<Real> cE(z.cE), r1(z.r1), r2(z.r2), lE(z.lE), cI(z.cI), s1(z.s1), s2(z.s2), lI(z.lI);
    Vector<Real> x(z.x);

    // Backtracking loop
    while (a.p >= EPSILON)
    {
      // Set trial point
      updatePoint();
      evalFunctions();

      // Check for function evaluation error
      if (z.err == 0)
      {
        // Set remaining trial values
        evalSlacks();
        evalMerit();

        // Check for nonlinear fraction-to-boundary violation
        Integer ftb{0};
        const Real tmp_min{std::min(p.ls_frac, z.mu)};
        if (i.nE > 0) {
          ftb += (z.r1 < tmp_min*r1).count() + (z.r2 < tmp_min*r2).count();
        }
        if (i.nI > 0) {
          ftb += (z.s1 < tmp_min*s1).count() + (z.s2 < tmp_min*s2).count();
        }

        // Check Armijo condition
        if (ftb == 0 && z.phi - phi <= -p.ls_thresh*a.p*std::max(d.qtred, 0.0))
        {
          // Reset variables and return
          setPrimals( x, r1, r2, s1, s2, lE, lI, z.f, z.cE, z.cI, z.phi);
          return;
        }
        else
        {
          // Reset variables
          setPrimals( x, r1, r2, s1, s2, lE, lI, f, cE, cI, phi);
        }
      }
      else
      {
        // Reset variables
        setPrimals( x, r1, r2, s1, s2, lE, lI, f, cE, cI, phi);
      }

      // Reduce step length
      a.p *= p.ls_factor;
    }
  }

  /**
   * \brief Fraction-to-boundary rule for primal and dual variables.
   *
   * Computes the maximum allowable primal and dual step fractions that keep interior-point slacks
   * and multipliers within their prescribed bounds. The computed values are stored in the \p a.p0
   * (initial primal fraction), \p a.p (primal steplength) and \p a.d (dual steplength).
   * \tparam Real Floating-point type used by the algorithm.
   * \param[in] a Acceptance object updated with computed fraction values.
   * \param[in] p Algorithm parameters.
   * \param[in] i Problem input structure.
   * \param[in] z Current iterate (used for slack/multiplier values).
   * \param[in] d Candidate search direction (may be scaled by the caller).
   */
  template <typename Real>
  inline
  void
  Solver<Real>::fractionToBoundary() {

    // Create alias for easier access
    Parameter<Real>  & p{this->m_parameter};
    Input<Real>      & i{this->m_input};
    Iterate<Real>    & z{this->m_iterate};
    Direction<Real>  & d{this->m_direction};
    Acceptance<Real> & a{this->m_acceptance};
    
    // Initialize primal fraction-to-boundary
    a.p0 = 1.0;

    // Update primal fraction-to-boundary for constraint slacks
    Real frac{std::min(p.ls_frac, z.mu)};
    if (i.nE > 0) {
      const Indices idx_r1(find(d.r1 < 0.0)), idx_r2(find(d.r2 < 0.0));
      Real min_1{std::numeric_limits<Real>::infinity()},
           min_2{std::numeric_limits<Real>::infinity()};
      if (idx_r1.size() > 0) {min_1 = (((frac - 1.0)*z.r1(idx_r1))/d.r1(idx_r1)).minCoeff();}
      if (idx_r2.size() > 0) {min_2 = (((frac - 1.0)*z.r2(idx_r2))/d.r2(idx_r2)).minCoeff();}
      a.p0 = std::min(a.p0, std::min(min_1, min_2));
    }
    if (i.nI > 0) {
      const Indices idx_s1(find(d.s1 < 0.0)), idx_s2(find(d.s2 < 0.0));
      Real min_1{std::numeric_limits<Real>::infinity()},
           min_2{std::numeric_limits<Real>::infinity()};
      if (idx_s1.size() > 0) {min_1 = (((frac - 1.0)*z.s1(idx_s1))/d.s1(idx_s1)).minCoeff();}
      if (idx_s2.size() > 0) {min_2 = (((frac - 1.0)*z.s2(idx_s2))/d.s2(idx_s2)).minCoeff();}
      a.p0 = std::min(a.p0, std::min(min_1, min_2));
    }

    // Initialize primal step length
    a.p = a.p0;

    // Initialize dual fraction-to-boundary
    a.d = 1.0;

    // Update dual fraction-to-boundary for constraint multipliers
    if (i.nE > 0) {
      const Indices idx_l(find(d.lE < 0.0)), idx_g(find(d.lE > 0));
      Real min_1{std::numeric_limits<Real>::infinity()},
           min_2{std::numeric_limits<Real>::infinity()};
      if (idx_l.size() > 0) {min_1 = (((frac - 1.0)*(1.0 + z.lE(idx_l)))/d.lE(idx_l)).minCoeff();}
      if (idx_g.size() > 0) {min_2 = (((1.0 - frac)*(1.0 - z.lE(idx_g)))/d.lE(idx_g)).minCoeff();}
      a.d = std::min(a.d, std::min(min_1, min_2));
    }
    if (i.nI > 0) {
      const Indices idx_l(find(d.lI < 0.0)), idx_g(find(d.lI > 0));
      Real min_1{std::numeric_limits<Real>::infinity()},
           min_2{std::numeric_limits<Real>::infinity()};
      if (idx_l.size() > 0) {min_1 = (((frac - 1.0)*       z.lI(idx_l)) /d.lI(idx_l)).minCoeff();}
      if (idx_g.size() > 0) {min_2 = (((1.0 - frac)*(1.0 - z.lI(idx_g)))/d.lI(idx_g)).minCoeff();}
      a.d = std::min(a.d, std::min(min_1, min_2));
    }
  }

  /**
   * \brief Full-step check for trial penalty parameters.
   *
   * Attempts full steps for a sequence of trial penalty parameters, verifying the Armijo condition
   * on the merit. If a full step satisfies the acceptance condition the routine returns 1 and the
   * iterate is left at the accepted trial point; otherwise the penalty parameter is decreased and
   * the process repeats.
   * \tparam Real Floating-point type used by the algorithm.
   * \param[in] a Acceptance object with current fractions.
   * \param[in] p Algorithm parameters.
   * \param[in] i Problem input structure.
   * \param[in] c Counters used during evaluations.
   * \param[in] z Current iterate (modified while testing trial points).
   * \param[in] d Candidate search direction.
   * \param[in] problem Problem interface used to evaluate objective/constraints.
   * \return 1 if a full step was accepted, 0 otherwise.
   */
  template <typename Real>
  inline
  Integer
  Solver<Real>::fullStepCheck() {

    // Create alias for easier access
    Parameter<Real>  & p{this->m_parameter};
    Input<Real>      & i{this->m_input};
    Iterate<Real>    & z{this->m_iterate};
    Direction<Real>  & d{this->m_direction};
    Acceptance<Real> & a{this->m_acceptance};

    // Set current and last penalty parameters
    Real rho{z.rho}, rho_temp{z.rho_};

    // Loop through last penalty parameters
    while (rho < rho_temp)
    {
      // Set penalty parameter
      setRho(rho_temp);

      // Evaluate merit
      evalMerit();

      // Store current values
      Real phi{z.phi}, f{z.f};
      Array<Real> cE(z.cE), r1(z.r1), r2(z.r2), lE(z.lE), cI(z.cI), s1(z.s1), s2(z.s2), lI(z.lI);
      Vector<Real> x(z.x);

      // Set trial point
      updatePoint();
      evalFunctions();

      // Check for function evaluation error
      if (z.err == 0)
      {
        // Set remaining trial values
        evalSlacks();
        evalMerit();

        // Check for nonlinear fraction-to-boundary violation
        Integer ftb{0};
        const Real tmp_min{std::min(p.ls_frac, z.mu)};
        if (i.nE > 0) {
          ftb += (z.r1 < tmp_min*r1).count() + (z.r2 < tmp_min*r2).count();
        }
        if (i.nI > 0) {
          ftb += (z.s1 < tmp_min*s1).count() + (z.s2 < tmp_min*s2).count();
        }

        // Check Armijo condition
        if (ftb == 0 && z.phi - phi <= -p.ls_thresh*a.p*std::max(d.qtred, 0.0))
        {
          // Reset variables, set boolean, and return
          setPrimals( x, r1, r2, s1, s2, lE, lI, z.f, z.cE, z.cI, z.phi);
          return 1;
        }
        else
        {
          // Reset variables
          setPrimals( x, r1, r2, s1, s2, lE, lI, f, cE, cI, phi);
        }
      }
      else
      {
        // Reset variables
        setPrimals( x, r1, r2, s1, s2, lE, lI, f, cE, cI, phi);
      }

      // Decrease rho
      rho_temp *= p.rho_factor;
    }

    // Set rho
    setRho(rho);

    // Evaluate merit
    evalMerit();

    return 0;
  }

  /**
   * \brief Line-search driver.
   *
   * Runs the fraction-to-boundary rule, full-step check, second-order correction and backtracking
   * as needed to find an acceptable trial step. The acceptance object is updated with the final
   * chosen step.
   * \tparam Real Floating-point type used by the algorithm.
   * \param[in] a Acceptance object updated by the routine.
   * \param[in] p Algorithm parameters.
   * \param[in] i Problem input structure.
   * \param[in] c Counters used during evaluations.
   * \param[in] z Current iterate (modified temporarily while testing points).
   * \param[in] d Candidate search direction (modified as trial directions are formed).
   * \param[in] problem Problem interface used to evaluate objective/constraints.
   */
  template <typename Real>
  inline
  void
  Solver<Real>::lineSearch() {

    // Create alias for easier access
    Acceptance<Real> & a{this->m_acceptance};

    // Check fraction-to-boundary rule
    fractionToBoundary();

    // Check for full step for trial penalty parameters
    Integer b{ fullStepCheck() };
    // Run second-order correction
    a.s = false;
    if (b == 0) {
      b = secondOrderCorrection();
      if (b == 2) {a.s = true;}
    }

    // Run backtracking line search
    if (b == 0) backtracking();
  }

  /**
   * \brief Second-order correction (SOC).
   *
   * Attempts a second-order correction to improve feasibility when the trial step fails the Armijo
   * condition (or other checks). The routine evaluates a corrected Newton step and tests acceptance
   * similarly to the main line-search path.
   * \tparam Real Floating-point type used by the algorithm.
   * \param[in] a Acceptance object containing current step fractions (updated on success).
   * \param[in] p Algorithm parameters.
   * \param[in] i Problem input structure.
   * \param[in] c Counters used during evaluations.
   * \param[in] z Current iterate (modified while testing trial points).
   * \param[in] d Candidate search direction (modified to hold trial direction).
   * \param[in] problem Problem interface used to evaluate objective/constraints.
   * \return 0 if SOC not accepted, 1 if first acceptance criterion met, 2 if SOC accepted.
   */
  template <typename Real>
  inline
  Integer
  Solver<Real>::secondOrderCorrection() {

    // Create alias for easier access
    Parameter<Real>  & p{this->m_parameter};
    Input<Real>      & i{this->m_input};
    Iterate<Real>    & z{this->m_iterate};
    Direction<Real>  & d{this->m_direction};
    Acceptance<Real> & a{this->m_acceptance};

    // Store current iterate values
    Real f{z.f}, phi{z.phi}, v{z.v};
    Array<Real> cE(z.cE), r1(z.r1), r2(z.r2), lE(z.lE), cI(z.cI), s1(z.s1), s2(z.s2), lI(z.lI);
    Vector<Real> x(z.x);

    // Set trial point
    updatePoint();
    evalFunctions();

    // Check for function evaluation error
    if (z.err == 0)
    {
      // Set remaining trial values
      evalSlacks();
      evalMerit();

      // Check for nonlinear fraction-to-boundary violation
      Integer ftb{0};
      const Real tmp_min{std::min(p.ls_frac, z.mu)};
      if (i.nE > 0) {
        ftb += (z.r1 < tmp_min*r1).count() + (z.r2 < tmp_min*r2).count();
      }
      if (i.nI > 0) {
        ftb += (z.s1 < tmp_min*s1).count() + (z.s2 < tmp_min*s2).count();
      }

      // Check Armijo condition
      if (ftb == 0 && z.phi - phi <= -p.ls_thresh*a.p*std::max(d.qtred, 0.0))
      {
        // Reset variables, set flag, and return
        setPrimals( x, r1, r2, s1, s2, lE, lI, z.f, z.cE, z.cI, z.phi);
        return 1;
      }
      else if ( evalViolation(z.cE, z.cI) < v )
      {
        // Reset variables and return
        setPrimals( x, r1, r2, s1, s2, lE, lI, f, cE, cI, phi);
        return 0;
      }
      else
      {
        // Reset variables (but leave constraint values for second-order correction)
        setPrimals( x, r1, r2, s1, s2, lE, lI, f, z.cE, z.cI, phi);
      }
    }
    else
    {
      // Reset variables and return
      setPrimals( x, r1, r2, s1, s2, lE, lI, f, cE, cI, phi);
      return 0;
    }

    // Recompute slacks for second order correction
    evalSlacks();

    // Evaluate trial primal-dual right-hand side vector
    evalNewtonRhs();

    // Store current direction values
    Vector<Real> dx(d.x), dr1(d.r1), dr2(d.r2), dlE(d.lE), ds1(d.s1), ds2(d.s2), dlI(d.lI);
    Real dx_norm{d.x_norm}, dl_norm{d.l_norm};

    // Evaluate search direction
    evalNewtonStep();

    // Set trial direction
    setDirection(
      a.p*dx+d.x,
      a.p*dr1+d.r1.matrix(),
      a.p*dr2+d.r2.matrix(),
      a.p*ds1+d.s1.matrix(),
      a.p*ds2+d.s2.matrix(),
      a.d*dlE+d.lE.matrix(),
      a.d*dlI+d.lI.matrix(),
      (a.p*dx+d.x).norm(),
      std::sqrt( (a.d*dlE+d.lE.matrix()).matrix().squaredNorm() +
                 (a.d*dlI+d.lI.matrix()).squaredNorm() )
    );

    // Set trial point
    updatePoint();
    evalFunctions();

    // Check for function evaluation error
    if (z.err == 0)
    {
      // Set remaining trial values
      evalSlacks();
      evalMerit();

      // Check for nonlinear fraction-to-boundary violation
      Integer ftb{0};
      const Real tmp_min{std::min(p.ls_frac, z.mu)};
      if (i.nE > 0) {
        ftb += (z.r1 < tmp_min*r1).count() + (z.r2 < tmp_min*r2).count();
      }
      if (i.nI > 0) {
        ftb += (z.s1 < tmp_min*s1).count() + (z.s2 < tmp_min*s2).count();
      }

      // Check Armijo condition
      if (ftb == 0 && z.phi - phi <= -p.ls_thresh*a.p*std::max(d.qtred, 0.0))
      {
        // Reset variables, set flag, and return
        setPrimals( x, r1, r2, s1, s2, lE, lI, z.f, z.cE, z.cI, z.phi);
        return 2;
      }
      else
      {
        // Reset variables
        setPrimals( x, r1, r2, s1, s2, lE, lI, f, cE, cI, phi);
      }
    }
    else
    {
      // Reset variables
      setPrimals( x, r1, r2, s1, s2, lE, lI, f, cE, cI, phi);
    }

    // Reset direction
    setDirection(dx, dr1, dr2, ds1, ds2, dlE, dlI, dx_norm, dl_norm);

    // Reduce step length
    a.p *= p.ls_factor;

    return 0;
  }

} // namespace Pipal

#endif /* INCLUDE_PIPAL_ACCEPTANCE_HH */
