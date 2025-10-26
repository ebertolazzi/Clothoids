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

#ifndef INCLUDE_PIPAL_OPTIMIZE_HH
#define INCLUDE_PIPAL_OPTIMIZE_HH

namespace Pipal {

  template <typename Real>
  inline
  bool
  Solver<Real>::optimize( Vector<Real> const & x_guess, Vector<Real> & x_sol ) {
    
    #define CMD "Pipal::Solver::optimize(...): "

    // Create alias for easier access
    Counter          & c{this->m_counter};
    Input<Real>      & i{this->m_input};
    Iterate<Real>    & z{this->m_iterate};
    Direction<Real>  & d{this->m_direction};
    Acceptance<Real> & a{this->m_acceptance};

    // Check that the problem is set
    PIPAL_ASSERT(
      this->m_problem.get() != nullptr,
      CMD "problem not set, use 'problem(...)' method to set it"
    );

    // Get variable bounds
    Vector<Real> bl, bu;
    PIPAL_ASSERT(
      this->m_problem->primal_lower_bounds(bl),
      CMD "error in evaluating lower bounds on primal variables"
    );
    PIPAL_ASSERT(
      this->m_problem->primal_upper_bounds(bu),
      CMD "error in evaluating upper bounds on primal variables"
    );

    // Get constraint bounds
    Vector<Real> cl, cu;
    PIPAL_ASSERT(
      this->m_problem->constraints_lower_bounds(cl),
      CMD "error in evaluating lower bounds on constraints"
    );
    PIPAL_ASSERT(
      this->m_problem->constraints_upper_bounds(cu),
      CMD "error in evaluating upper bounds on constraints"
    );

    // Reset counters
    resetCounter();

    // Fill input structure
    buildInput(this->m_problem->name(), x_guess, bl, bu, cl, cu);
    buildIterate();
    resetDirection(d);

    // Print header and break line
    if (this->m_verbose) {this->m_output.printHeader(i, z); this->m_output.printBreak(c);}

    // Iterations loop
    while (!checkTermination()) {

      // Print iterate
      if (this->m_verbose) {this->m_output.printIterate(c, z);}

      // Evaluate the step
      evalStep();

      // Print direction
      if (this->m_verbose) {this->m_output.printDirection(z, d);}

      lineSearch();

      // Print accepted
      if (this->m_verbose) {this->m_output.printAcceptance(a);}

      updateIterate();

      // Increment iteration counter
      incrementIterationCount();

      // Print break
      if (this->m_verbose) {this->m_output.printBreak(c);}
    }
    // Print footer and terminate
    if (this->m_verbose) {this->m_output.printFooter(c, z, checkTermination() );}

    // Get solution in original variables
    evalXOriginal(x_sol);

    // Return success if is finite
    return x_sol.allFinite();

    #undef CMD
  }

} // namespace Solver

#endif /* INCLUDE_PIPAL_OPTIMIZE_HXX */
