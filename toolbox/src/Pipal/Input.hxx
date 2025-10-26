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

#ifndef INCLUDE_PIPAL_INPUT_HXX
#define INCLUDE_PIPAL_INPUT_HXX

namespace Pipal
{

  // Constructor
  template <typename Real>
  inline
  void
  Solver<Real>::buildInput(
    std::string  const & name,
    Vector<Real> const & x0,
    Vector<Real> const & bl,
    Vector<Real> const & bu,
    Vector<Real> const & cl,
    Vector<Real> const & cu
  ) {

    // Create alias for easier access
    Parameter<Real> & p{this->m_parameter};
    Input<Real>     & i{this->m_input};

    //#define CMD "Pipal::resetInput(...): "

    // Set problem identity
    i.name = name;

    // Set number of original formulation variables
    i.n0 = static_cast<Integer>(x0.size());

    // Find indices sets
    static constexpr Real EPSILON{std::numeric_limits<Real>::epsilon()};
    Mask const cond_bl(bl.array() <= -p.rhs_bnd);
    Mask const cond_bu(bu.array() >= p.rhs_bnd);
    Mask const cond_cl(cl.array() <= -p.rhs_bnd);
    Mask const cond_cu(cu.array() >= p.rhs_bnd);
    Mask const cond_beq((bl.array() - bu.array()).abs() <= EPSILON);
    Mask const cond_ceq((cl.array() - cu.array()).abs() <= EPSILON);

    // Now the index sets
    i.I1 = find(cond_bl && cond_bu);
    i.I2 = find(cond_beq);
    i.I3 = find(!cond_bl && cond_bu);
    i.I4 = find(cond_bl && !cond_bu);
    i.I5 = find(!cond_bl && !cond_bu && !cond_beq);
    i.I6 = find(cond_ceq);
    i.I7 = find(!cond_cl && cond_cu);
    i.I8 = find(cond_cl && !cond_cu);
    i.I9 = find(!cond_cl && !cond_cu && !cond_ceq);

    // Set right-hand side values
    i.b2 = bl(i.I2);
    i.l3 = bl(i.I3);
    i.u4 = bu(i.I4);
    i.l5 = bl(i.I5);
    i.u5 = bu(i.I5);
    i.b6 = cl(i.I6);
    i.l7 = cl(i.I7);
    i.u8 = cu(i.I8);
    i.l9 = cl(i.I9);
    i.u9 = cu(i.I9);

    // Set sizes of indices sets
    i.n1 = static_cast<Integer>(i.I1.size());
    i.n2 = static_cast<Integer>(i.I2.size());
    i.n3 = static_cast<Integer>(i.I3.size());
    i.n4 = static_cast<Integer>(i.I4.size());
    i.n5 = static_cast<Integer>(i.I5.size());
    i.n6 = static_cast<Integer>(i.I6.size());
    i.n7 = static_cast<Integer>(i.I7.size());
    i.n8 = static_cast<Integer>(i.I8.size());
    i.n9 = static_cast<Integer>(i.I9.size());

    // Initialize number of invalid bounds
    i.vi = 0;

    // Count invalid bounds
    if (i.n2 > 0) {
      i.vi += (i.b2.array() <= -p.rhs_bnd).count();
      i.vi += (i.b2.array() >= p.rhs_bnd).count();
    }
    if (i.n3 > 0) {
      i.vi += (i.l3.array() >= p.rhs_bnd).count();
    }
    if (i.n4 > 0) {
      i.vi += (i.u4.array() <= -p.rhs_bnd).count();
    }
    if (i.n5 > 0) {
      i.vi += (i.l5.array() >= p.rhs_bnd).count();
      i.vi += (i.u5.array() <= -p.rhs_bnd).count();
      i.vi += (i.l5.array() > i.u5.array()).count();
    }
    if (i.n6 > 0) {
      i.vi += (i.b6.array() <= -p.rhs_bnd).count();
      i.vi += (i.b6.array() >= p.rhs_bnd).count();
    }
    if (i.n7 > 0) {
      i.vi += (i.l7.array() >= p.rhs_bnd).count();
    }
    if (i.n8 > 0) {
      i.vi += (i.u8.array() <= -p.rhs_bnd).count();
    }
    if (i.n9 > 0) {
      i.vi += (i.l9.array() >= p.rhs_bnd).count();
      i.vi += (i.u9.array() <= -p.rhs_bnd).count();
      i.vi += (i.l9.array() > i.u9.array()).count();
    }

    // Set number of variables and constraints
    i.nV = i.n1 + i.n3 + i.n4 + i.n5;
    i.nI = i.n3 + i.n4 + 2*i.n5 + i.n7 + i.n8 + 2*i.n9;
    i.nE = i.n6;

    // Set size of primal-dual matrix
    i.nA = i.nV + 3*i.nE + 3*i.nI;

    // Set initial point
    i.x0.resize(i.nV);
    i.x0 << x0(i.I1), x0(i.I3), x0(i.I4), x0(i.I5);

    //#undef CMD
  }

} // namespace Pipal

#endif /* INCLUDE_PIPAL_INPUT_HH */
