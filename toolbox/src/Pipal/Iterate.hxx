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

#ifndef INCLUDE_PIPAL_ITERATE_HXX
#define INCLUDE_PIPAL_ITERATE_HXX

namespace Pipal {

  /**
   * \brief Initialize an Iterate object for a given problem/input.
   * \tparam Real Floating-point type used by the algorithm.
   * \param[out] z Iterate object to initialize.
   * \param[in] p Algorithm parameters (may be read).
   * \param[in] i Problem input structure.
   * \param[in] c Counters used to account evaluations.
   * \param[in] problem Problem interface used for objective/constraints.
   */
  template <typename Real>
  inline
  void
  Solver<Real>::buildIterate() {

    // Create alias for easier access
    Parameter<Real> & p{this->m_parameter};
    Input<Real>     & i{this->m_input};
    Iterate<Real>   & z{this->m_iterate};

    // Safe initialization of all members
    z.f = 0;
    z.fu = 0;
    z.g.setZero(i.nV);
    z.r1.setZero(i.nE);
    z.r2.setZero(i.nE);
    z.cE.setZero(i.nE);
    z.JE.resize(i.nE, i.nV);
    z.s1.setZero(i.nI);
    z.s2.setZero(i.nI);
    z.cI.setZero(i.nI);
    z.JI.resize(i.nI, i.nV);
    z.H.resize(i.nV, i.nV);
    z.v = 0;
    z.vu = 0;
    z.phi = 0;
    z.Annz = 0;
    z.shift = 0;
    z.b.resize(i.nA);
    z.kkt.setZero(3);
    z.kkt_.setConstant(p.opt_err_mem, std::numeric_limits<Real>::infinity());
    z.fs = 1.0;
    z.cEs.setOnes(i.nE);
    z.cEu.setZero(i.nE);
    z.cIs.setOnes(i.nI);
    z.cIu.setZero(i.nI);
    z.A.resize(i.nA, i.nA);
    z.shift22 = 0;
    z.cut_ = false;

    // Initialize point
    z.x     = i.x0;
    z.rho   = p.rho_init;
    z.mu    = p.mu_init;
    z.lE.setZero(i.nE);
    z.lI.setConstant(i.nI, 0.5);
    z.err   = 0;
    evalScalings();
    evalFunctions();
    evalGradients();
    evalDependent();
    z.v0    = 1.0;
    evalInfeasibility(z);
    z.v0    = z.v;
    evalInfeasibility(z);
    z.v_    = z.v;
    evalHessian();
    z.Hnnz  = static_cast<Integer>(z.H.nonZeros());
    z.JEnnz = static_cast<Integer>(z.JE.nonZeros());
    z.JInnz = static_cast<Integer>(z.JI.nonZeros());
    initNewtonMatrix();
    evalNewtonMatrix();
  }

  /**
   * \brief Check termination criteria for the solver.
   *
   * Evaluates optimality/feasibility/iteration-count/bounds error conditions and returns the
   * termination code used by the solver.
   * \tparam Real Floating-point type used by the algorithm.
   * \return An integer termination code (0 = continue, 1..5 = specific exits).
   */
  template <typename Real>
  inline
  Integer
  Solver<Real>::checkTermination() const {

    // Create alias for easier access
    Parameter<Real> const & p{this->m_parameter};
    Counter         const & c{this->m_counter};
    Input<Real>     const & i{this->m_input};
    Iterate<Real>   const & z{this->m_iterate};

    // Update termination based on optimality error of nonlinear optimization problem
    if (z.kkt(1) <= p.opt_err_tol && z.v <= p.opt_err_tol) return 1;

    // Update termination based on optimality error of feasibility problem
    if (z.kkt(0) <= p.opt_err_tol && z.v > p.opt_err_tol) return 2;

    // Update termination based on iteration count
    if (c.k >= p.iter_max) return 3;

    // Update termination based on invalid bounds
    if (i.vi > 0) return 4;

    // Update termination based on function evaluation error
    if (z.err > 0) return 5;

    return 0;
  }

  /**
   * \brief Evaluate objective and constraint functions at the current iterate.
   * \tparam Real Floating-point type used by the algorithm.
   * \param[in] z Current iterate.
   * \param[in] i Problem input structure.
   * \param[in] c Counters used to account evaluations.
   * \param[in] problem Problem interface used for objective/constraints.
   */
  template <typename Real>
  inline
  void
  Solver<Real>::evalFunctions() {
  
    // Create alias for easier access
    Input<Real>   & i{this->m_input};
    Iterate<Real> & z{this->m_iterate};

    // Evaluate x in original space
    Vector<Real> x_orig;
    evalXOriginal(x_orig);

    // Initialize/Reset evaluation flag
    z.err = 0;

    // Increment function evaluation counter
    incrementFunctionCount();

    // Try AMPL functions evaluation
    Vector<Real> c_orig;
    try
    {
      // Evaluate AMPL functions
      m_problem->objective(x_orig, z.f);
      m_problem->constraints(x_orig, c_orig);
    }
    catch (...)
    {
      // Set evaluation flag, default values, and return
      z.err = 1;
      z.f   = std::numeric_limits<Real>::quiet_NaN();
      z.cE.setConstant(i.nE, std::numeric_limits<Real>::quiet_NaN());
      z.cI.setConstant(i.nI, std::numeric_limits<Real>::quiet_NaN());
      z.fu  = std::numeric_limits<Real>::quiet_NaN();
      z.cEu.setConstant(i.nE, std::numeric_limits<Real>::quiet_NaN());
      z.cIu.setConstant(i.nI, std::numeric_limits<Real>::quiet_NaN());
      return;
    }

    // Set equality constraint values
    if (i.nE > 0) {z.cE = c_orig(i.I6) - i.b6;}

    // Initialize inequality constraint values
    if (i.nI > 0) {z.cI.setZero(i.nI);}

    // Set inequality constraint values
    if (i.n3 > 0) {
      z.cI.head(i.n3) = i.l3 - z.x.segment(i.n1, i.n3);
    }
    if (i.n4 > 0) {
      z.cI.segment(i.n3, i.n4) = -i.u4 + z.x.segment(i.n1+i.n3, i.n4);
    }
    if (i.n5 > 0) {
      z.cI.segment(i.n3+i.n4,      i.n5) =  i.l5 - z.x.segment(i.n1+i.n3+i.n4, i.n5);
      z.cI.segment(i.n3+i.n4+i.n5, i.n5) = -i.u5 + z.x.segment(i.n1+i.n3+i.n4, i.n5);
    }
    if (i.n7 > 0) {
      z.cI.segment(i.n3+i.n4+i.n5+i.n5, i.n7) = i.l7 - c_orig(i.I7);
    }
    if (i.n8 > 0) {
      z.cI.segment(i.n3+i.n4+i.n5+i.n5+i.n7, i.n8) = -i.u8 + c_orig(i.I8);
    }
    if (i.n9 > 0) {
      z.cI.segment(i.n3+i.n4+i.n5+i.n5+i.n7+i.n8,      i.n9) =  i.l9 - c_orig(i.I9);
      z.cI.segment(i.n3+i.n4+i.n5+i.n5+i.n7+i.n8+i.n9, i.n9) = -i.u9 + c_orig(i.I9);
    }

    // Store unscaled quantities
    z.fu = z.f;
    if (i.nE > 0) {z.cEu = z.cE;}
    if (i.nI > 0) {z.cIu = z.cI;}

    // Scale quantities
    z.f *= z.fs;
    if (i.nE > 0) {z.cE = (z.cEs*z.cE).matrix();}
    if (i.nI > 0) {z.cI = (z.cIs*z.cI).matrix();}

  }


  /**
   * \brief Insert a dense block into a sparse matrix at the specified offsets.
   * \tparam Real Floating-point type used by the algorithm.
   * \tparam CheckZero Whether to skip zero entries when inserting.
   * \param[out] mat_sparse Sparse matrix where to insert the block.
   * \param[in] mat_dense Dense matrix block to insert.
   * \param[in] row_offset Row offset in the sparse matrix.
   * \param[in] col_offset Column offset in the sparse matrix.
   */
  template <typename Real, bool CheckZero = true>
  inline
  void
  insert_block(
    SparseMatrix<Real> & mat_sparse,
    Matrix<Real> const & mat_dense,
    Integer      const   row_offset,
    Integer      const   col_offset
  ) {
    const Integer cols{static_cast<Integer>(mat_dense.cols())}, rows{static_cast<Integer>(mat_dense.rows())};
    for (Integer r{0}; r < rows; ++r) {
      for (Integer c{0}; c < cols; ++c) {
        //if constexpr (CheckZero) {
        //  if (mat_dense(r, c) != static_cast<Real>(0.0)) {
        //    mat_sparse.coeffRef(r + row_offset, c + col_offset) = mat_dense(r, c);
        //  }
        //} else {
        mat_sparse.coeffRef(r + row_offset, c + col_offset) = mat_dense(r, c);
        //}
      }
    }
  }

  /**
   * \brief Evaluate objective gradient and constraint Jacobian.
   *
   * Calls the problem-provided gradient/jacobian functions in the original space, maps the results
   * into the internal compressed representations and applies scaling. Increments gradient counters
   * and handles exceptions.
   * \tparam Real Floating-point type used by the algorithm.
   * \param[in] z Current iterate.
   * \param[in] i Problem input structure.
   * \param[in] c Counters used to account evaluations.
   * \param[in] problem Problem interface used for objective/constraints.
   */
  template <typename Real>
  inline
  void
  Solver<Real>::evalGradients() {

    // Create alias for easier access
    Input<Real>   & i{this->m_input};
    Iterate<Real> & z{this->m_iterate};

    // Evaluate x in original space
    Vector<Real> x_orig;
    evalXOriginal(x_orig);

    // Initialize/Reset evaluation flag
    z.err = 0;

    // Increment gradient evaluation counter
    incrementGradientCount();

    // Try AMPL gradients evaluation
    Vector<Real> g_orig;
    Matrix<Real> J_orig;
    try
    {
      // Evaluate AMPL gradients
      m_problem->objective_gradient(x_orig, g_orig);
      SparseMatrix<Real> J_orig_sparse;
      m_problem->constraints_jacobian(x_orig, J_orig_sparse);
      J_orig = J_orig_sparse;
    }
    catch (...)
    {
      // Set evaluation flag, default values, and return
      z.err = 1;
      return;
    }

    // Set objective gradient
    z.g << g_orig(i.I1), g_orig(i.I3), g_orig(i.I4), g_orig(i.I5);

    // Set equality constraint Jacobian
    if (i.nE > 0) {
      Integer col_offset{0};
      insert_block<Real>(z.JE, J_orig(i.I6, i.I1), 0, col_offset);
      col_offset += i.I1.size();
      insert_block<Real>(z.JE, J_orig(i.I6, i.I3), 0, col_offset);
      col_offset += i.I3.size();
      insert_block<Real>(z.JE, J_orig(i.I6, i.I4), 0, col_offset);
      col_offset += i.I4.size();
      insert_block<Real>(z.JE, J_orig(i.I6, i.I5), 0, col_offset);
    }

    // Initialize inequality constraint Jacobian
    //if (i.nI > 0) {z.JI.resize(i.nI, i.nV);}

    // Set inequality constraint Jacobian
    if (i.n3 > 0) {
      for (Integer k{0}; k < i.n3; ++k)  {z.JI.coeffRef(k, i.n1+k) = -1.0;}
    }
    if (i.n4 > 0) {
      Integer tmp{i.n1+i.n3};
      for (Integer k{0}; k < i.n4; ++k) {z.JI.coeffRef(i.n3+k, tmp+k) = 1.0;}
    }
    if (i.n5 > 0) {
      Integer tmp_row{i.n3+i.n4}, tmp_col{i.n1+i.n3+i.n4};
      for (Integer k{0}; k < i.n5; ++k) {
        z.JI.coeffRef(tmp_row+k,      tmp_col+k) = -1.0;
        z.JI.coeffRef(tmp_row+i.n5+k, tmp_col+k) = 1.0;
      }
    }
    if (i.n7 > 0) {
      Integer row_offset{i.n3+i.n4+i.n5+i.n5}, col_offset{0};
      insert_block<Real>(z.JI, -J_orig(i.I7, i.I1), row_offset, col_offset);
      col_offset += i.I1.size();
      insert_block<Real>(z.JI, J_orig(i.I7, i.I3), row_offset, col_offset);
      col_offset += i.I3.size();
      insert_block<Real>(z.JI, J_orig(i.I7, i.I4), row_offset, col_offset);
      col_offset += i.I4.size();
      insert_block<Real>(z.JI, J_orig(i.I7, i.I5), row_offset, col_offset);
    }
    if (i.n8 > 0) {
      Integer row_offset{i.n3+i.n4+i.n5+i.n5+i.n7}, col_offset{0};
      insert_block<Real>(z.JI, J_orig(i.I8, i.I1), row_offset, col_offset);
      col_offset += i.I1.size();
      insert_block<Real>(z.JI, J_orig(i.I8, i.I3), row_offset, col_offset);
      col_offset += i.I3.size();
      insert_block<Real>(z.JI, J_orig(i.I8, i.I4), row_offset, col_offset);
      col_offset += i.I4.size();
      insert_block<Real>(z.JI, J_orig(i.I8, i.I5), row_offset, col_offset);
    }
    if (i.n9 > 0) {
      Integer row_offset{i.n3+i.n4+i.n5+i.n5+i.n7+i.n8}, col_offset{0};
      insert_block<Real>(z.JI, -J_orig(i.I9, i.I1), row_offset, col_offset);
      col_offset += i.I1.size();
      insert_block<Real>(z.JI, -J_orig(i.I9, i.I3), row_offset, col_offset);
      col_offset += i.I3.size();
      insert_block<Real>(z.JI, -J_orig(i.I9, i.I4), row_offset, col_offset);
      col_offset += i.I4.size();
      insert_block<Real>(z.JI, -J_orig(i.I9, i.I5), row_offset, col_offset);
      row_offset += i.n9; // next row block
      col_offset = 0;
      insert_block<Real>(z.JI, J_orig(i.I9, i.I1), row_offset, col_offset);
      col_offset += i.I1.size();
      insert_block<Real>(z.JI, J_orig(i.I9, i.I3), row_offset, col_offset);
      col_offset += i.I3.size();
      insert_block<Real>(z.JI, J_orig(i.I9, i.I4), row_offset, col_offset);
      col_offset += i.I4.size();
      insert_block<Real>(z.JI, J_orig(i.I9, i.I5), row_offset, col_offset);
    }

    // Scale objective gradient
    z.g *= z.fs;

    // Scale constraint Jacobians
    if (i.nE > 0) {
      for (Integer k{0}; k < z.JE.outerSize(); ++k) {
        for (typename SparseMatrix<Real>::InnerIterator it(z.JE, k); it; ++it) {
          it.valueRef() *= z.cEs[it.row()];
        }
      }
    }
    if (i.nI > 0) {
      for (Integer k{0}; k < z.JI.outerSize(); ++k) {
        for (typename SparseMatrix<Real>::InnerIterator it(z.JI, k); it; ++it) {
          it.valueRef() *= z.cIs[it.row()];
        }
      }
    }
  }

  /**
   * \brief Evaluate the Hessian of the Lagrangian and assemble internal H.
   *
   * Calls the problem's Hessian-of-Lagrangian routine (with multipliers and primal variables in the
   * original space), maps blocks into the internal sparse Hessian and rescales according to
   * penalization and scaling factors.
   * \tparam Real Floating-point type used by the algorithm.
   * \param[in] z Current iterate.
   * \param[in] i Problem input structure.
   * \param[in] c Counters used to account evaluations.
   * \param[in] problem Problem interface used for objective/constraints.
   */
  template <typename Real>
  inline
  void
  Solver<Real>::evalHessian() {
    // Create alias for easier access
    Input<Real>   & i{this->m_input};
    Iterate<Real> & z{this->m_iterate};

    // Evaluate lambda in original space
    Vector<Real> l_orig, x_orig;
    evalLambdaOriginal(l_orig);
    evalXOriginal(x_orig);

    // Initialize/Reset evaluation flag
    z.err = 0;

    // Increment Hessian evaluation counter
    incrementHessianCount();

    // Try AMPL Hessian evaluation
    Matrix<Real> H_orig;
    try
    {
      // Evaluate H_orig
      SparseMatrix<Real> H_orig_sparse;
      m_problem->lagrangian_hessian(x_orig, l_orig, H_orig_sparse);
      H_orig = H_orig_sparse;
    }
    catch (...)
    {
      // Set evaluation flag, default values, and return
      z.err = 1;
      return;
    }

    // Set Hessian of the Lagrangian
    Integer row_offset{0}, col_offset{0};
    insert_block<Real>(z.H, H_orig(i.I1, i.I1), row_offset, col_offset);
    col_offset += i.I1.size();
    insert_block<Real>(z.H, H_orig(i.I1, i.I3), row_offset, col_offset);
    col_offset += i.I3.size();
    insert_block<Real>(z.H, H_orig(i.I1, i.I4), row_offset, col_offset);
    col_offset += i.I4.size();
    insert_block<Real>(z.H, H_orig(i.I1, i.I5), row_offset, col_offset);
    row_offset += i.I1.size();
    col_offset = 0; // next row block
    insert_block<Real>(z.H, H_orig(i.I3, i.I1), row_offset, col_offset);
    col_offset += i.I1.size();
    insert_block<Real>(z.H, H_orig(i.I3, i.I3), row_offset, col_offset);
    col_offset += i.I3.size();
    insert_block<Real>(z.H, H_orig(i.I3, i.I4), row_offset, col_offset);
    col_offset += i.I4.size();
    insert_block<Real>(z.H, H_orig(i.I3, i.I5), row_offset, col_offset);
    row_offset += i.I3.size();
    col_offset = 0; // next row block
    insert_block<Real>(z.H, H_orig(i.I4, i.I1), row_offset, col_offset);
    col_offset += i.I1.size();
    insert_block<Real>(z.H, H_orig(i.I4, i.I3), row_offset, col_offset);
    col_offset += i.I3.size();
    insert_block<Real>(z.H, H_orig(i.I4, i.I4), row_offset, col_offset);
    col_offset += i.I4.size();
    insert_block<Real>(z.H, H_orig(i.I4, i.I5), row_offset, col_offset);
    row_offset += i.I4.size();
    col_offset = 0; // next row block
    insert_block<Real>(z.H, H_orig(i.I5, i.I1), row_offset, col_offset);
    col_offset += i.I1.size();
    insert_block<Real>(z.H, H_orig(i.I5, i.I3), row_offset, col_offset);
    col_offset += i.I3.size();
    insert_block<Real>(z.H, H_orig(i.I5, i.I4), row_offset, col_offset);
    col_offset += i.I4.size();
    insert_block<Real>(z.H, H_orig(i.I5, i.I5), row_offset, col_offset);
    
    // workaround se mancano elementi sulla diagonale
    for ( Integer i{0}; i < z.H.rows(); ++i ) (void) z.H.coeffRef(i, i);  // accede o crea H(i,i)

    // Rescale H
    z.H *= z.rho*z.fs;
  }

  /**
   * \brief Compute the infinity-norm of the KKT optimality vector.
   * \tparam Real Floating-point type used by the algorithm.
   * \param[in] z Current iterate.
   * \param[in] i Problem input structure.
   * \param[in] rho Current penalty parameter.
   * \param[in] mu Current interior-point parameter.
   * \return The infinity-norm of the KKT optimality vector.
   */
  template <typename Real>
  inline
  Real
  Solver<Real>::evalKKTError(
    Real const rho,
    Real const mu
  ) {
    // Create alias for easier access
    Input<Real>   & i{this->m_input};
    Iterate<Real> & z{this->m_iterate};

    // Initialize optimality vector
    Vector<Real> kkt(i.nV+2*i.nE+2*i.nI);
    kkt.setZero();

    // Set gradient of penalty objective
    kkt.head(i.nV) = rho*z.g;

    // Set gradient of Lagrangian for constraints
    if (i.nE > 0) {kkt.head(i.nV) += (z.lE.matrix().transpose()*z.JE).transpose();}
    if (i.nI > 0) {kkt.head(i.nV) += (z.lI.matrix().transpose()*z.JI).transpose();}

    // Set complementarity for constraint slacks
    if (i.nE > 0) {
      kkt.segment(i.nV, 2*i.nE) << z.r1*(1.0 + z.lE) - mu, z.r2*(1.0 - z.lE) - mu;
    }
    if (i.nI > 0) {
      kkt.segment(i.nV+2*i.nE, 2*i.nI) << z.s1*z.lI - mu, z.s2*(1.0 - z.lI) - mu;
    }

    // Scale complementarity
    if (rho > 0) {
      kkt *= (1.0 / std::max(1.0, (rho*z.g).template lpNorm<Eigen::Infinity>()));
    }

    // Evaluate optimality error
    return kkt.template lpNorm<Eigen::Infinity>();
  }

  /**
   * \brief Compute the three KKT error measures used by the solver.
   * \tparam Real Floating-point type used by the algorithm.
   * \param[in] z Current iterate.
   * \param[in] i Problem input structure.
   */
  template <typename Real>
  inline
  void
  Solver<Real>::evalKKTErrors() {
    // Create alias for easier access
    Iterate<Real> & z{this->m_iterate};

    // Loop to compute optimality errors
    z.kkt(0) = evalKKTError(0.0,   0.0);
    z.kkt(1) = evalKKTError(z.rho, 0.0);
    z.kkt(2) = evalKKTError(z.rho, z.mu);
  }

  /**
   * \brief Reconstruct multipliers in the original variable/constraint space.
   * \tparam Real Floating-point type used by the algorithm.
   * \param[in] z Current iterate.
   * \param[in] i Problem input structure.
   * \param[out] l Vector to store multipliers in original space.
   */
  template <typename Real>
  inline
  void
  Solver<Real>::evalLambdaOriginal( Vector<Real> & l ) const {
    // Create alias for easier access
    Input<Real>   const & i{this->m_input};
    Iterate<Real> const & z{this->m_iterate};
    
    // Initialize multipliers in original space
    l.setZero(i.nE+i.n7+i.n8+i.n9);

    // Scale equality constraint multipliers
    if (i.nE > 0) {l(i.I6) = z.lE*(z.cEs/(z.rho*z.fs));}

    // Scale inequality constraint multipliers
    Array<Real> lI;
    if (i.n7+i.n8+i.n9 > 0) {lI = z.lI*(z.cIs/(z.rho*z.fs));}

    // Set inequality constraint multipliers in original space
    if (i.n7 > 0) {
      l(i.I7) = -lI.segment(i.n3+i.n4+i.n5+i.n5, i.n7);
    }
    if (i.n8 > 0) {
      l(i.I8) = lI.segment(i.n3+i.n4+i.n5+i.n5+i.n7, i.n8);
    }
    if (i.n9 > 0) {
      l(i.I9) = lI.segment(i.n3+i.n4+i.n5+i.n5+i.n7+i.n8+i.n9, i.n9)
                -lI.segment(i.n3+i.n4+i.n5+i.n5+i.n7+i.n8, i.n9);
    }
  }


  /**
   * \brief Compute the merit function value for the current iterate.
   * \tparam Real Floating-point type used by the algorithm.
   * \param[in] z Current iterate.
   * \param[in] i Problem input structure.
   */
  template <typename Real>
  inline
  void
  Solver<Real>::evalMerit() {

    // Create alias for easier access
    Input<Real>   & i{this->m_input};
    Iterate<Real> & z{this->m_iterate};

    // Initialize merit for objective
    z.phi = z.rho*z.f;

    // Update merit for slacks
    if (i.nE > 0) {
      Array<Real> r_all(2*i.nE); r_all << z.r1, z.r2;
      z.phi -= z.mu * r_all.log().sum() - r_all.sum();
    }
    if (i.nI > 0) {
      Array<Real> s_all(2*i.nI); s_all << z.s1, z.s2;
      z.phi -= z.mu * s_all.log().sum() - z.s2.sum();
    }
  }

  /**
   * \brief Assemble and (attempt to) factorize the Newton system matrix.
   * \tparam Real Floating-point type used by the algorithm.
   * \param[in] z Current iterate.
   * \param[in] p Algorithm parameters.
   * \param[in] i Problem input structure.
   * \param[in] c Counters used to account evaluations.
   */
  template <typename Real>
  inline
  void
  Solver<Real>::evalNewtonMatrix() {

    // Create alias for easier access
    Parameter<Real> & p{this->m_parameter};
    Input<Real>     & i{this->m_input};
    Iterate<Real>   & z{this->m_iterate};

    // Check for equality constraints
    if (i.nE > 0)
    {
      // Set diagonal terms
      for (Integer j{0}; j < i.nE; ++j) {
        z.A.coeffRef(i.nV+j,      i.nV+j)      = (1.0 + z.lE(j))/z.r1(j);
        z.A.coeffRef(i.nV+i.nE+j, i.nV+i.nE+j) = (1.0 - z.lE(j))/z.r2(j);
      }

      // Set constraint Jacobian
      insert_block<Real>(z.A, z.JE, i.nV+2*i.nE+2*i.nI, 0);
    }

    // Check for inequality constraints
    if (i.nI > 0)
    {
      // Set diagonal terms
      const Integer offset{i.nV+2*i.nE};
      for (Integer j{0}; j < i.nI; ++j) {
        z.A.coeffRef(offset+j,      offset+j)      = z.lI(j)/z.s1(j);
        z.A.coeffRef(offset+i.nI+j, offset+i.nI+j) = (1.0 - z.lI(j))/z.s2(j);
      }

      // Set constraint Jacobian
      insert_block<Real>(z.A, z.JI, i.nV+3*i.nE+2*i.nI, 0);
    }

    // Set minimum potential shift
    Real min_shift{std::max(p.shift_min, p.shift_factor1*z.shift)};

    // Initialize Hessian modification
    if (z.cut_ == true) {z.shift = std::min(p.shift_max, min_shift/p.shift_factor2);}
    else {z.shift = 0;}
    // Initialize inertia correction loop
    bool done{false};
    z.shift22 = 0;

    // Loop until inertia is correct
    while (!done && z.shift < p.shift_max)
    {
      // Set Hessian of Lagrangian
      insert_block<Real>(z.A, z.H+z.shift*Matrix<Real>::Identity(i.nV, i.nV), 0, 0);

      // Set diagonal terms
      {
        const Integer offset{i.nV+2*i.nE+2*i.nI};
        for (Integer j{0}; j < i.nE; ++j) {
          z.A.coeffRef(offset+j, offset+j) = -z.shift22;
        }
      }
      {
        const Integer offset{i.nV+3*i.nE+2*i.nI};
        // Set diagonal terms
        for (Integer j{0}; j < i.nI; ++j) {
          z.A.coeffRef(offset+j, offset+j) = -z.shift22;
        }
      }

      // Set number of nonzeros in (upper triangle of) Newton matrix
      z.Annz = static_cast<Integer>(z.A.nonZeros());

      // Factor primal-dual matrix
      z.ldlt.compute(z.A);

      // Approximate number of negative pivots (inertia)
      Integer neig{static_cast<Integer>((z.ldlt.vectorD().array() < 0.0).count())};

      // Increment factorization counter
      incrementFactorizationCount();

      // Set number of nonnegative eigenvalues
      Integer peig{i.nA - neig};

      // Check inertia
      if (peig < i.nV+2*i.nE+2*i.nI) {z.shift = std::max(min_shift, z.shift/p.shift_factor2);}
      else if (neig < i.nE+i.nI && z.shift22 == 0) {z.shift22 = p.shift_min;}
      else {done = true;}
    }

    // Update Hessian
    z.H.diagonal().array() += z.shift;

    // Compress Newton matrix
    z.A.makeCompressed();
  }

  /**
   * \brief Build the right-hand side vector for the Newton system.
   * \tparam Real Floating-point type used by the algorithm.
   * \param[in] z Current iterate.
   * \param[in] i Problem input structure.
   */
  template <typename Real>
  inline
  void
  Solver<Real>::evalNewtonRhs() {
    // Create alias for easier access
    Input<Real>   & i{this->m_input};
    Iterate<Real> & z{this->m_iterate};

    // Initialize right-hand side vector
    z.b.setZero();

    // Set gradient of objective
    for (Integer k{0}; k < i.nV; ++k) {
      for (Integer k{0}; k < i.nV; ++k) {z.b.coeffRef(k) = z.rho*z.g(k);}

      // Set gradient of Lagrangian for constraints
      if (i.nE > 0) {
        Vector<Real> tmp((z.lE.matrix().transpose()*z.JE).transpose());
        for (Integer k{0}; k < tmp.size(); ++k) {z.b.coeffRef(k) += tmp[k];}
      }
      if (i.nI > 0) {
        Vector<Real> tmp((z.lI.matrix().transpose()*z.JI).transpose());
        for (Integer k{0}; k < tmp.size(); ++k) {z.b.coeffRef(k) += tmp[k];}
      }
    }

    // Set complementarity for constraint slacks
    if (i.nE > 0) {
      // Compute element-wise complementarity terms with safe element-wise division
      Vector<Real> tmp(2*i.nE);
      tmp << 1.0 + z.lE - z.mu * z.r1.cwiseInverse(), 1.0 - z.lE - z.mu * z.r2.cwiseInverse();
      for (Integer k{0}; k < 2*i.nE; ++k) {z.b.insert(i.nV + k) = tmp[k];}
    }
    if (i.nI > 0) {
      // Compute element-wise complementarity terms with safe element-wise division
      Vector<Real> tmp(2*i.nI);
      tmp << z.lI - z.mu * z.s1.cwiseInverse(), 1.0 - z.lI - z.mu * z.s2.cwiseInverse();
      for (Integer k{0}; k < 2*i.nI; ++k) {z.b.insert(i.nV + 2*i.nE + k) = tmp[k];}
    }

    // Set penalty-interior-point constraint values
    if (i.nE > 0) {
      const Vector<Real> tmp(z.cE + z.r1 - z.r2);
      const Integer offset{i.nV+2*i.nE+2*i.nI};
      for (Integer k{0}; k < i.nE; ++k) {z.b.coeffRef(offset+k) = tmp[k];}
    }
    if (i.nI > 0) {
      const Vector<Real> tmp(z.cI + z.s1 - z.s2);
      const Integer offset{i.nV+3*i.nE+2*i.nI};
      for (Integer k{0}; k < i.nI; ++k) {z.b.coeffRef(offset+k) = tmp[k];}
    }
  }

  /**
   * \brief Evaluate scaling multipliers for objective and constraints.
   * \tparam Real Floating-point type used by the algorithm.
   * \param[in] z Current iterate.
   * \param[in] p Algorithm parameters.
   * \param[in] i Problem input structure.
   * \param[in] c Counters used to account evaluations.
   * \param[in] problem Problem interface used for objective/constraints.
   */
  template <typename Real>
  inline
  void
  Solver<Real>::evalScalings() {

    // Create alias for easier access
    Parameter<Real> & p{this->m_parameter};
    Input<Real>     & i{this->m_input};
    Iterate<Real>   & z{this->m_iterate};

    // Initialize scalings
    z.fs = 1;
    z.cEs.setOnes(i.nE);
    z.cIs.setOnes(i.nI);

    // Evaluate gradients
    evalGradients();

    // Scale down objective if norm of gradient is too large
    z.fs = p.grad_max / std::max(z.g.template lpNorm<Eigen::Infinity>(), p.grad_max);

    // Loop through equality constraints
    for (Integer j{0}; j < i.nE; ++j)
    {
      // Scale down equality constraint j if norm of gradient is too large
      Real row_inf_norm{0.0};
      for (Integer c{0}; c < z.JE.outerSize(); ++c) {
        for (typename SparseMatrix<Real>::InnerIterator it(z.JE, c); it; ++it) {
          if (it.row() == j) {row_inf_norm = std::max(row_inf_norm, std::abs(it.value()));}
        }
      }
      z.cEs(j) = p.grad_max / std::max(row_inf_norm, p.grad_max);
    }

    // Loop through inequality constraints
    for (Integer j{0}; j < i.nI; ++j)
    {
      // Scale down inequality constraint j if norm of gradient is too large
      Real row_inf_norm{0.0};
      for (Integer c{0}; c < z.JI.outerSize(); ++c) {
        for (typename SparseMatrix<Real>::InnerIterator it(z.JI, c); it; ++it) {
          if (it.row() == j) {row_inf_norm = std::max(row_inf_norm, std::abs(it.value()));}
        }
      }
      z.cIs(j) = p.grad_max / std::max(row_inf_norm, p.grad_max);
    }
  }

  /**
   * \brief Compute internal slack variables from current iterate.
   * \tparam Real Floating-point type used by the algorithm.
   * \param[in] z Current iterate.
   * \param[in] p Algorithm parameters.
   * \param[in] i Problem input structure.
   */
  template <typename Real>
  inline
  void
  Solver<Real>::evalSlacks() {
    // Create alias for easier access
    Parameter<Real> & p{this->m_parameter};
    Input<Real>     & i{this->m_input};
    Iterate<Real>   & z{this->m_iterate};
 
    // Check for equality constraints
    if (i.nE > 0)
    {
      // Set slacks
      z.r1 = 0.5*((z.mu - z.cE) + (z.cE.square() + z.mu*z.mu).sqrt());
      z.r2 = 0.5*((z.mu + z.cE) + (z.cE.square() + z.mu*z.mu).sqrt());

      // Adjust for numerical error
      z.r1.cwiseMax(p.slack_min);
      z.r2.cwiseMax(p.slack_min);
    }

    // Check for inequality constraints
    if (i.nI > 0)
    {
      // Set slacks
      z.s1 = 0.5*((2.0*z.mu - z.cI) + (z.cI.square() + 4.0*z.mu*z.mu).sqrt());
      z.s2 = 0.5*((2.0*z.mu + z.cI) + (z.cI.square() + 4.0*z.mu*z.mu).sqrt());

      // Adjust for numerical error
      z.s1.cwiseMax(p.slack_min);
      z.s2.cwiseMax(p.slack_min);
    }
  }

  /**
   * \brief Reconstruct the full primal vector in the original variable ordering.
   * \tparam Real Floating-point type used by the algorithm.
   * \param[in] z Current iterate.
   * \param[in] i Problem input structure.
   * \param[out] x Vector to store primal variables in original space.
   */
  template <typename Real>
  inline
  void
  Solver<Real>::evalXOriginal( Vector<Real> & x ) {

    // Create alias for easier access
    Input<Real>      & i{this->m_input};
    Iterate<Real>    & z{this->m_iterate};

    // Initialize x in original space
    x.setZero(i.n0);

    // Evaluate x in original space
    x(i.I1) = z.x.head(i.n1);
    x(i.I2) = i.b2;
    x(i.I3) = z.x.segment(i.n1, i.n3);
    x(i.I4) = z.x.segment(i.n1+i.n3, i.n4);
    x(i.I5) = z.x.segment(i.n1+i.n3+i.n4, i.n5);
  }

  /**
   * \brief Reserve and initialize the internal sparse Newton matrix structure.
   * \tparam Real Floating-point type used by the algorithm.
   * \param[in] z Current iterate.
   * \param[in] i Problem input structure.
   */
  template <typename Real>
  inline
  void
  Solver<Real>::initNewtonMatrix() {

    // Create alias for easier access
    Input<Real>   & i{this->m_input};
    Iterate<Real> & z{this->m_iterate};

    // Allocate memory
    z.A.reserve(z.Hnnz + 5*i.nE + 5*i.nI + z.JEnnz + z.JInnz);

    // Initialize interior-point Hessians
    {
      Integer diag;
      for (Integer k{0}; k < 2*i.nE; ++k) {
        diag = i.nV+k;
        z.A.coeffRef(diag, diag) = 1.0;
      }

      for (Integer k{0}; k < 2*i.nI; ++k) {
        diag = i.nV+2*i.nE+k;
        z.A.coeffRef(diag, diag) = 1.0;
      }
    }

    // Check for constraints
    if (i.nE > 0)
    {
      // Initialize constraint Jacobian
      Integer row, col;
      for (Integer k{0}; k < i.nE; ++k) {
        row = i.nV+2*i.nE+2*i.nI+k, col = i.nV+k;
        z.A.coeffRef(row, col) = 1.0;
        z.A.coeffRef(row, col+i.nE) = -1.0;
      }
    }

    // Check for inequality constraints
    if (i.nI > 0)
    {
      // Initialize constraint Jacobian
      Integer row, col;
      for (Integer k{0}; k < i.nI; ++k) {
        row = i.nV+3*i.nE+2*i.nI+k, col = i.nV+2*i.nE+k;
        z.A.coeffRef(row, col) = 1.0;
        z.A.coeffRef(row, col+i.nI) = -1.0;
      }
    }
  }

  /**
   * \brief Set primal/dual blocks and associated quantities on an iterate.
   * \tparam Real Floating-point type used by the algorithm.
   * \param[in] z Current iterate.
   * \param[in] i Problem input structure.
   * \param[in] x Primal variable vector.
   * \param[in] r1 Equality constraint slack variables (lower).
   * \param[in] r2 Equality constraint slack variables (upper).
   * \param[in] s1 Inequality constraint slack variables (lower).
   * \param[in] s2 Inequality constraint slack variables (upper).
   * \param[in] lE Equality constraint multipliers.
   * \param[in] lI Inequality constraint multipliers.
   * \param[in] f Objective function value.
   * \param[in] cE Equality constraint values.
   * \param[in] cI Inequality constraint values.
   * \param[in] phi Merit function value.
   */
  template <typename Real>
  inline
  void
  Solver<Real>::setPrimals(
    Vector<Real> const & x,
    Array<Real>  const & r1,
    Array<Real>  const & r2,
    Array<Real>  const & s1,
    Array<Real>  const & s2,
    Array<Real>  const & lE,
    Array<Real>  const & lI,
    Real         const   f,
    Array<Real>  const & cE,
    Array<Real>  const & cI,
    Real         const   phi
  ) {
    // Create alias for easier access
    Input<Real>   & i{this->m_input};
    Iterate<Real> & z{this->m_iterate};

    // Set primal variables
    z.x = x; z.f = f;
    if (i.nE > 0) {z.cE = cE; z.r1 = r1; z.r2 = r2; z.lE = lE;}
    if (i.nI > 0) {z.cI = cI; z.s1 = s1; z.s2 = s2; z.lI = lI;}
    z.phi = phi;
  }

  /**
   * \brief Update the iterate after a trial step is accepted.
   *
   * Applies the accepted step to the iterate, recomputes infeasibility and gradients and updates
   * parameter memory used for adaptive strategies.
   * \tparam Real Floating-point type used by the algorithm.
   * \param[in] z Current iterate.
   * \param[in] p Algorithm parameters.
   * \param[in] i Problem input structure.
   * \param[in] c Counters used to account evaluations.
   * \param[in] d Computed search direction.
   * \param[in] a Step acceptance information.
   * \param[in] problem Problem interface used for objective/constraints.
   */
  template <typename Real>
  inline
  void
  Solver<Real>::updateIterate() {

    // Create alias for easier access
    Parameter<Real>  & p{this->m_parameter};
    Iterate<Real>    & z{this->m_iterate};
    Acceptance<Real> & a{this->m_acceptance};

    // Update last quantities
    z.v_   = z.v;
    z.cut_ = (a.p < a.p0);

    // Update iterate quantities
    updatePoint();
    evalInfeasibility(z);
    evalGradients();
    evalDependent();

    // Update last KKT errors
    //z.kkt_.resize(p.opt_err_mem);
    z.kkt_ << z.kkt(1), z.kkt_.head(p.opt_err_mem-1);
  }

  /**
   * \brief Update penalty and interior-point parameters based on KKT errors.
   *
   * Adjusts rho and mu using the solver's adaptive rules to drive optimality and feasibility towards
   * the desired tolerances.
   * \tparam Real Floating-point type used by the algorithm.
   * \param[in] z Current iterate.
   * \param[in] p Algorithm parameters.
   * \param[in] i Problem input structure.
   */
  template <typename Real>
  inline
  void
  Solver<Real>::updateParameters() {

    Parameter<Real> & p{this->m_parameter};
    Iterate<Real>   & z{this->m_iterate};

    // Check for interior-point parameter update based on optimality error
    while (z.mu > p.mu_min && z.kkt(2) <= std::max(z.mu, p.opt_err_tol-z.mu))
    {
      // Restrict interior-point parameter increase
      setMuMaxExpZero();

      // Update interior-point parameter
      if (z.mu > p.mu_min)
      {
        // Decrease interior-point
        z.mu = std::max(p.mu_min, std::min(p.mu_factor*z.mu, std::pow(z.mu, p.mu_factor_exp)));

        // Evaluate penalty and interior-point parameter dependent quantities
        evalDependent();
      }
    }

    // Check for penalty parameter update based on optimality error
    if ((z.kkt(1) <= p.opt_err_tol && z.v > p.opt_err_tol) ||
          z.v > std::max({1.0, z.v_, p.infeas_max}))
    {
      // Update penalty parameter
      if (z.rho > p.rho_min)
      {
        // Decrease penalty parameter
        z.rho = std::max(p.rho_min, p.rho_factor*z.rho);

        // Evaluate penalty and interior-point parameter dependent quantities
        evalDependent();
      }
    }
  }

  /**
   * \brief Apply a step to the primal and dual variables of the iterate.
   * \tparam Real Floating-point type used by the algorithm.
   * \param[in] z Current iterate.
   * \param[in] i Problem input structure.
   * \param[in] d Computed search direction.
   * \param[in] a Step acceptance information.
   */
  template <typename Real>
  inline
  void
  Solver<Real>::updatePoint() {
    // Create alias for easier access
    Input<Real>      & i{this->m_input};
    Iterate<Real>    & z{this->m_iterate};
    Direction<Real>  & d{this->m_direction};
    Acceptance<Real> & a{this->m_acceptance};

    // Update primal and dual variables
    z.x += a.p*d.x ;
    if (i.nE > 0) {z.r1 += a.p*d.r1; z.r2 += a.p*d.r2;}
    if (i.nI > 0) {z.s1 += a.p*d.s1; z.s2 += a.p*d.s2;}
    if (i.nE > 0) {z.lE += a.d*d.lE;}
    if (i.nI > 0) {z.lI += a.d*d.lI;}
  }

  /**
   * \brief Compute the 1-norm feasibility violation from equality/inequality values.
   * \tparam Real Floating-point type used by the algorithm.
   * \param[in] i Problem input structure.
   * \param[in] cE Equality constraint values.
   * \param[in] cI Inequality constraint values.
   * \return The 1-norm feasibility violation.
   */
  template <typename Real>
  inline
  Real
  Solver<Real>::evalViolation(
    Array<Real> const & cE,
    Array<Real> const & cI
  ) const {
    // Create alias for easier access
    Input<Real> const & i{this->m_input};
  
    // Initialize violation vector
    Vector<Real> vec;

    // Update vector for constraint values
    if (i.nE > 0) {vec = cE;}
    if (i.nI > 0) {
      Vector<Real> cIpos(cI.cwiseMax(0.0));
      if (vec.size() > 0) {
        Vector<Real> tmp(vec.size() + cIpos.size());
        tmp << vec, cIpos;
        vec = std::move(tmp);
      } else {
        vec = cIpos;
      }
    }

    // Evaluate vector norm
    return vec.template lpNorm<1>();
  }

} // namespace Pipal

#endif /* INCLUDE_PIPAL_ITERATE_HH */
