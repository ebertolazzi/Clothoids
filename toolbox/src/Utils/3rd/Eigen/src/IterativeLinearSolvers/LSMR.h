// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-FileCopyrightText: The Eigen Authors
// SPDX-License-Identifier: MPL-2.0

#ifndef EIGEN_LSMR_H
#define EIGEN_LSMR_H

// LSMR is an iterative algorithm for least-squares problems  min ||A x - b||,
// based on the Golub-Kahan bidiagonalization process.  It is analytically
// equivalent to applying MINRES to the normal equation  A^T A x = A^T b, so the
// quantity ||A^T r_k|| decreases monotonically (where r_k = b - A x_k).  In
// practice ||r_k|| also decreases monotonically, which makes LSMR safer than
// LSQR to stop early.  With a damping parameter lambda > 0 it instead minimizes
// || (A; lambda I) x - (b; 0) ||, i.e. it solves the regularized (Tikhonov)
// least-squares problem.
//
// This implementation follows the published algorithm:
//
//   D. C.-L. Fong and M. A. Saunders, "LSMR: An Iterative Algorithm for Sparse
//   Least-Squares Problems", SIAM J. Sci. Comput. 33(5):2950-2971, 2011.
//   https://web.stanford.edu/group/SOL/software/lsmr/
//
// The scalar recurrences, the cheap estimates of ||r||, ||A^T r||, ||A|| and
// cond(A), and the stopping rules reproduce the reference Fortran 90
// implementation by the same authors (distributed under the BSD / Common Public
// License at the URL above).  No source code from that implementation is
// reproduced here; only the mathematical algorithm is implemented, in Eigen's
// own idiom.

// IWYU pragma: private
#include "./InternalHeaderCheck.h"

namespace Eigen {

namespace internal {

/** \internal Low-level LSMR algorithm
 * \param mat       The matrix A.
 * \param rhs       The right hand side vector b.
 * \param x         On input the initial guess x0 (usually zero), on output the
 *                  computed solution.
 * \param precond   A preconditioner.  It is applied as a self-adjoint right
 *                  preconditioner: LSMR is run on the better-conditioned system
 *                  \f$ A M^{-1} z = b \f$ (where \c precond.solve() applies
 *                  \f$ M^{-1} \f$) and the solution is recovered as
 *                  \f$ x = M^{-1} z \f$.  With IdentityPreconditioner this
 *                  reduces exactly to the reference algorithm.
 * \param iters     On input the maximum number of iterations, on output the
 *                  number of iterations performed.
 * \param tol_error On output, an estimate of the relative residual of the
 *                  normal equations \f$ ||A^T r|| / (||A||\,||r||) \f$.
 * \param atol      Stopping tolerance bounding the assumed relative error in
 *                  the entries of \a A.
 * \param btol      Stopping tolerance bounding the assumed relative error in
 *                  the entries of \a b.
 * \param lambda    The damping parameter \f$ \lambda \ge 0 \f$ (0 for an
 *                  unregularized problem).
 * \param conlim    An upper limit on cond(A); iterations stop if the estimate
 *                  exceeds it.  Pass 0 to disable (equivalent to 1/eps).
 *
 * \returns the LSMR stopping code \c istop:
 *   - 0: x = 0 (i.e. the supplied guess) is the exact solution; no iterations.
 *   - 1: A x = b is compatible and ||r|| is small enough, given atol and btol.
 *   - 2: a least-squares solution is good enough, given atol.
 *   - 3: the estimate of cond(A) exceeded conlim.
 *   - 4: A x = b is compatible and ||r|| is as small as machine precision allows.
 *   - 5: a least-squares solution is as good as machine precision allows.
 *   - 6: cond(A) seems too large for this machine.
 *   - 7: the iteration limit was reached.
 */
template <typename MatrixType, typename Rhs, typename Dest, typename Preconditioner>
EIGEN_DONT_INLINE Index lsmr(const MatrixType& mat, const Rhs& rhs, Dest& x, const Preconditioner& precond,
                             Index& iters, typename Dest::RealScalar& tol_error, const typename Dest::RealScalar& atol,
                             const typename Dest::RealScalar& btol, const typename Dest::RealScalar& lambda,
                             const typename Dest::RealScalar& conlim) {
  using numext::abs;
  using numext::sqrt;
  typedef typename Dest::RealScalar RealScalar;
  typedef typename Dest::Scalar Scalar;
  typedef Matrix<Scalar, Dynamic, 1> VectorType;

  const RealScalar zero(0);
  const RealScalar one(1);

  const Index n = mat.cols();
  const Index maxIters = iters;

  // n-vectors needed before the early-return below. u is the only m-vector.
  VectorType v(n), Atu(n);

  // Set up the first vectors u and v for the bidiagonalization. These satisfy
  // beta*u = b - A*x0  and  alpha*v = M^{-1} A^T u.
  VectorType u = rhs - mat * x;  // working residual r0 of the initial guess
  RealScalar alpha = zero;
  RealScalar beta = u.stableNorm();
  if (beta > zero) {
    u /= beta;
    Atu.noalias() = mat.adjoint() * u;
    v = precond.solve(Atu);
    alpha = v.stableNorm();
    if (alpha > zero) v /= alpha;
  }

  iters = 0;
  // If b - A*x0 = 0 or A^T(b - A*x0) = 0 then the current x already solves the
  // (least-squares) problem: no correction is needed.
  if (alpha * beta == zero) {
    tol_error = zero;
    return 0;  // istop = 0
  }

  // dx accumulates the correction (in z-space when a preconditioner is used).
  // The remaining n-vectors are only needed once we start iterating.
  VectorType dx = VectorType::Zero(n);
  VectorType h = v;
  VectorType hbar = VectorType::Zero(n);
  VectorType t(n);

  // Quantities driving the two plane rotations.
  RealScalar alphabar = alpha;
  RealScalar zetabar = alpha * beta;
  RealScalar rho = one;
  RealScalar rhobar = one;
  RealScalar cbar = one;
  RealScalar sbar = zero;

  // Quantities for the running estimate of ||r||.
  RealScalar betadd = beta;
  RealScalar betad = zero;
  RealScalar rhodold = one;
  RealScalar tautildeold = zero;
  RealScalar thetatilde = zero;
  RealScalar zeta = zero;
  RealScalar d = zero;

  // Quantities for the running estimates of ||A|| and cond(A).
  RealScalar normA2 = alpha * alpha;
  RealScalar maxrbar = zero;
  RealScalar minrbar = NumTraits<RealScalar>::highest();

  const RealScalar normb = beta;
  const RealScalar ctol = conlim > zero ? one / conlim : zero;
  RealScalar test2 = zero;  // recomputed each iteration; also read after the loop for tol_error

  Index istop = 0;
  Index itn = 0;
  while (istop == 0) {
    ++itn;

    // Perform the next step of the bidiagonalization to obtain the next
    // beta, u, alpha, v.  These satisfy
    //     beta*u = A*(M^{-1} v) - alpha*u,
    //     alpha*v = M^{-1} A^T u - beta*v.
    t = precond.solve(v);
    u *= -alpha;
    u.noalias() += mat * t;
    beta = u.stableNorm();
    if (beta > zero) {
      u /= beta;
      Atu.noalias() = mat.adjoint() * u;
      t = precond.solve(Atu);
      v = t - beta * v;
      alpha = v.stableNorm();
      if (alpha > zero) v /= alpha;
    }

    // Construct rotation Qhat_{k} that folds in the damping.
    const RealScalar alphahat = numext::hypot(alphabar, lambda);
    const RealScalar chat = alphabar / alphahat;
    const RealScalar shat = lambda / alphahat;

    // Use a plane rotation Q_{k} to turn B_{k} into R_{k}.
    const RealScalar rhoold = rho;
    rho = numext::hypot(alphahat, beta);
    const RealScalar c = alphahat / rho;
    const RealScalar s = beta / rho;
    const RealScalar thetanew = s * alpha;
    alphabar = c * alpha;

    // Use a plane rotation Qbar_{k} to turn R_{k}^T into Rbar_{k}.
    const RealScalar rhobarold = rhobar;
    const RealScalar zetaold = zeta;
    const RealScalar thetabar = sbar * rho;
    const RealScalar rhotemp = cbar * rho;
    rhobar = numext::hypot(cbar * rho, thetanew);
    cbar = cbar * rho / rhobar;
    sbar = thetanew / rhobar;
    zeta = cbar * zetabar;
    zetabar = -sbar * zetabar;

    // Update h, hbar and the (correction) solution dx.
    hbar = h - (thetabar * rho / (rhoold * rhobarold)) * hbar;
    dx += (zeta / (rho * rhobar)) * hbar;
    h = v - (thetanew / rho) * h;

    // Estimate ||r||.  Apply rotation Qhat_{k}, then Q_{k}, then Qtilde_{k-1}.
    const RealScalar betaacute = chat * betadd;
    const RealScalar betacheck = -shat * betadd;
    const RealScalar betahat = c * betaacute;
    betadd = -s * betaacute;

    const RealScalar thetatildeold = thetatilde;
    const RealScalar rhotildeold = numext::hypot(rhodold, thetabar);
    const RealScalar ctildeold = rhodold / rhotildeold;
    const RealScalar stildeold = thetabar / rhotildeold;
    thetatilde = stildeold * rhobar;
    rhodold = ctildeold * rhobar;
    betad = -stildeold * betad + ctildeold * betahat;

    tautildeold = (zetaold - thetatildeold * tautildeold) / rhotildeold;
    const RealScalar taud = (zeta - thetatilde * tautildeold) / rhodold;
    d += betacheck * betacheck;
    const RealScalar normr = sqrt(d + numext::abs2(betad - taud) + numext::abs2(betadd));

    // Estimate ||A||.
    normA2 += beta * beta;
    const RealScalar normA = sqrt(normA2);
    normA2 += alpha * alpha;

    // Estimate cond(A).
    maxrbar = numext::maxi(maxrbar, rhobarold);
    if (itn > 1) minrbar = numext::mini(minrbar, rhobarold);
    const RealScalar condA = numext::maxi(maxrbar, rhotemp) / numext::mini(minrbar, rhotemp);

    // Compute the norms needed for the stopping rules.
    const RealScalar normAr = abs(zetabar);
    const RealScalar normx = dx.stableNorm();

    const RealScalar test1 = normr / normb;
    test2 = (normA * normr > zero) ? normAr / (normA * normr) : zero;
    const RealScalar test3 = one / condA;
    const RealScalar t1 = test1 / (one + normA * normx / normb);
    const RealScalar rtol = btol + atol * normA * normx / normb;

    // The "1 + test <= 1" guards trigger near machine precision and make the
    // method behave as if atol = btol = eps and conlim = 1/eps even when the
    // user passed 0 for any of them.  The user tolerances are tested first so a
    // genuine convergence (istop 1/2/3) takes priority over the machine limits.
    // (istop 4 uses t1 rather than test1, matching the reference algorithm.)
    if (test1 <= rtol)
      istop = 1;
    else if (test2 <= atol)
      istop = 2;
    else if (test3 <= ctol)
      istop = 3;
    else if (one + t1 <= one)
      istop = 4;
    else if (one + test2 <= one)
      istop = 5;
    else if (one + test3 <= one)
      istop = 6;
    else if (itn >= maxIters)
      istop = 7;
  }

  // Recover the solution: x <- x0 + M^{-1} dx.
  t = precond.solve(dx);
  x += t;

  iters = itn;
  tol_error = test2;
  return istop;
}

}  // namespace internal

template <typename MatrixType_, typename Preconditioner_ = IdentityPreconditioner>
class LSMR;

namespace internal {

template <typename MatrixType_, typename Preconditioner_>
struct traits<LSMR<MatrixType_, Preconditioner_> > {
  typedef MatrixType_ MatrixType;
  typedef Preconditioner_ Preconditioner;
};

}  // namespace internal

/** \ingroup IterativeLinearSolvers_Module
 * \brief An LSMR solver for sparse (or dense) least-squares problems
 *
 * This class solves for the least-squares solution of \c A \c x = \c b using the
 * LSMR algorithm of Fong and Saunders.  LSMR is based on the Golub-Kahan
 * bidiagonalization and is analytically equivalent to MINRES applied to the
 * normal equation \f$ A^T A x = A^T b \f$, so the residual of the normal
 * equation \f$ ||A^T r|| \f$ decreases monotonically.  In exact arithmetic LSMR
 * returns the minimum-norm least-squares solution.  The matrix \c A can be
 * non-symmetric and rectangular; \c A and the vectors \c x and \c b can be
 * either dense or sparse.  LSMR only needs \c A through the products \f$ Av \f$
 * and \f$ A^T u \f$, so a matrix-free operator may be used as well.
 *
 * Unlike LeastSquaresConjugateGradient, which forms \f$ A^T A \f$ implicitly,
 * LSMR works on \c A directly through the bidiagonalization and is therefore
 * more robust on ill-conditioned problems.
 *
 * \tparam MatrixType_ the type of the matrix A, can be a dense or a sparse matrix.
 * \tparam Preconditioner_ the type of the preconditioner. Default is IdentityPreconditioner.
 *
 * \implsparsesolverconcept
 *
 * The maximum number of iterations and the tolerance can be controlled via the
 * setMaxIterations() and setTolerance() methods. The defaults are twice the
 * number of columns of the matrix for the maximum number of iterations and
 * NumTraits<Scalar>::epsilon() for the tolerance. setTolerance() sets both of
 * the algorithm's stopping tolerances \c atol (relative error assumed in \c A)
 * and \c btol (relative error assumed in \c b); they can also be set
 * independently via setToleranceA() and setToleranceB().
 *
 * The setDamping() method enables Tikhonov regularization: with a damping
 * \f$ \lambda > 0 \f$ the solver minimizes
 * \f$ ||Ax-b||^2 + \lambda^2 ||x||^2 \f$, for which a unique solution always
 * exists. The setConditionLimit() method can be used to stop the iterations as
 * soon as the estimated condition number of \c A exceeds a given bound.
 *
 * This class can be used like the other iterative solvers. Here is a typical
 * usage example:
   \code
   int m = 1000000, n = 10000;
   VectorXd x(n), b(m);
   SparseMatrix<double> A(m, n);
   // fill A and b
   LSMR<SparseMatrix<double> > lsmr;
   lsmr.compute(A);
   x = lsmr.solve(b);
   std::cout << "#iterations:     " << lsmr.iterations() << std::endl;
   std::cout << "estimated error: " << lsmr.error()      << std::endl;
   // update b, and solve again
   x = lsmr.solve(b);
   \endcode
 *
 * By default the iterations start with x=0 as an initial guess of the solution.
 * One can control the start using the solveWithGuess() method.
 *
 * If a non-default preconditioner is supplied it is applied as a self-adjoint
 * right preconditioner (LSMR is run on \f$ A M^{-1} z = b \f$ and the solution
 * recovered from \f$ M x = z \f$); it should therefore be self-adjoint, as
 * DiagonalPreconditioner and IdentityPreconditioner are.
 *
 * \sa class LeastSquaresConjugateGradient, class ConjugateGradient, SparseLU, SparseQR
 */
template <typename MatrixType_, typename Preconditioner_>
class LSMR : public IterativeSolverBase<LSMR<MatrixType_, Preconditioner_> > {
 protected:
  typedef IterativeSolverBase<LSMR> Base;
  using Base::m_error;
  using Base::m_info;
  using Base::m_isInitialized;
  using Base::m_iterations;
  using Base::matrix;

 public:
  typedef MatrixType_ MatrixType;
  typedef typename MatrixType::Scalar Scalar;
  typedef typename MatrixType::RealScalar RealScalar;
  typedef Preconditioner_ Preconditioner;

  /** Default constructor. */
  LSMR() : Base() {}

  /** Initialize the solver with matrix \a A for further \c Ax=b solving.
   *
   * This constructor is a shortcut for the default constructor followed
   * by a call to compute().
   *
   * \warning this class stores a reference to the matrix A as well as some
   * precomputed values that depend on it. Therefore, if \a A is changed
   * this class becomes invalid. Call compute() to update it with the new
   * matrix A, or modify a copy of A.
   */
  template <typename MatrixDerived>
  explicit LSMR(const EigenBase<MatrixDerived>& A) : Base(A.derived()) {}

  /** Sets the damping parameter \f$ \lambda \ge 0 \f$ for Tikhonov
   * regularization. With \a lambda > 0 the solver minimizes
   * \f$ ||Ax-b||^2 + \lambda^2 ||x||^2 \f$. The default is 0 (no damping).
   *
   * \note When a non-identity preconditioner \f$ M \f$ is used the damping
   * applies in the preconditioned variable, i.e. the solver minimizes
   * \f$ ||Ax-b||^2 + \lambda^2 ||Mx||^2 \f$. With the default
   * IdentityPreconditioner this is exactly \f$ ||Ax-b||^2 + \lambda^2 ||x||^2 \f$.
   */
  LSMR& setDamping(const RealScalar& lambda) {
    m_lambda = lambda;
    return *this;
  }

  /** \returns the damping parameter. \sa setDamping() */
  RealScalar damping() const { return m_lambda; }

  /** Sets an upper limit on the estimated condition number of \a A. The
   * iterations stop as soon as the estimate exceeds \a conlim. The default is
   * 0, which disables the limit (equivalent to 1/epsilon).
   */
  LSMR& setConditionLimit(const RealScalar& conlim) {
    m_conditionLimit = conlim;
    return *this;
  }

  /** \returns the condition-number limit. \sa setConditionLimit() */
  RealScalar conditionLimit() const { return m_conditionLimit; }

  /** Sets the stopping tolerance \c atol, which bounds the relative error
   * assumed in the entries of \a A. It drives the least-squares stopping rule
   * \f$ ||A^T r|| \le atol\,||A||\,||r|| \f$. If left unset (the default) it
   * falls back to tolerance(). \sa setToleranceB(), setTolerance() */
  LSMR& setToleranceA(const RealScalar& atol) {
    m_atol = atol;
    return *this;
  }

  /** \returns \c atol, or tolerance() if setToleranceA() has not been called.
   * \sa setToleranceA() */
  RealScalar toleranceA() const { return m_atol >= RealScalar(0) ? m_atol : Base::m_tolerance; }

  /** Sets the stopping tolerance \c btol, which bounds the relative error
   * assumed in the entries of \a b. It enters the compatible-system stopping
   * rule \f$ ||r|| \le btol\,||b|| + atol\,||A||\,||x|| \f$. If left unset (the
   * default) it falls back to tolerance(). \sa setToleranceA(), setTolerance() */
  LSMR& setToleranceB(const RealScalar& btol) {
    m_btol = btol;
    return *this;
  }

  /** \returns \c btol, or tolerance() if setToleranceB() has not been called.
   * \sa setToleranceB() */
  RealScalar toleranceB() const { return m_btol >= RealScalar(0) ? m_btol : Base::m_tolerance; }

  /** \internal */
  template <typename Rhs, typename Dest>
  void _solve_vector_with_guess_impl(const Rhs& b, Dest& x) const {
    m_iterations = Base::maxIterations();

    Index istop = internal::lsmr(matrix(), b, x, Base::m_preconditioner, m_iterations, m_error, toleranceA(),
                                 toleranceB(), m_lambda, m_conditionLimit);
    // istop in {0,1,2,4,5}: the (least-squares) solution was found, possibly
    // only to within machine precision (4,5). istop in {3,6,7}: stopped on the
    // condition-number limit or the iteration limit without meeting the
    // requested tolerance.
    m_info = (istop == 3 || istop == 6 || istop == 7) ? NoConvergence : Success;
  }

 protected:
  RealScalar m_lambda = RealScalar(0);
  RealScalar m_conditionLimit = RealScalar(0);
  // Negative means "unset": toleranceA()/toleranceB() then fall back to tolerance().
  RealScalar m_atol = RealScalar(-1);
  RealScalar m_btol = RealScalar(-1);
};

}  // end namespace Eigen

#endif  // EIGEN_LSMR_H
