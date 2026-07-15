// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2026 Rasmus Munk Larsen <rmlarsen@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0

#ifndef EIGEN_TRIDIAGONAL_INVERSE_ITERATION_H
#define EIGEN_TRIDIAGONAL_INVERSE_ITERATION_H

#include "./SelfAdjointEigenSolver.h"

// IWYU pragma: private
#include "./InternalHeaderCheck.h"

namespace Eigen {
namespace internal {

/** \internal
 *
 * Self-contained deterministic pseudo-random generator (SplitMix64) producing reproducible
 * values in the open interval (-1, 1). Used to seed the inverse-iteration start vectors.
 *
 * \c internal::random draws from the global \c rand() state, which is neither reproducible
 * across runs nor safe to call from several threads; a deterministic generator instead keeps
 * the eigenvector signs and the per-cluster reorthogonalization stable regardless of how the
 * columns are scheduled.
 */
struct inverse_iteration_rng {
  numext::uint64_t state;

  explicit inverse_iteration_rng(numext::uint64_t seed) : state(seed) {}

  template <typename RealScalar>
  RealScalar next() {
    // SplitMix64: a strong avalanche so that consecutive per-column seeds yield well-separated
    // (uncorrelated) start vectors -- consecutive LCG seeds would produce near-parallel starts that
    // collapse under the cluster reorthogonalization.
    state += 0x9E3779B97F4A7C15ULL;
    numext::uint64_t z = state;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    z = z ^ (z >> 31);
    // Top 32 bits give a uniform integer in [0, 2^32); map to [0, 1) and then to (-1, 1).
    const numext::uint32_t hi = numext::uint32_t(z >> 32);
    const RealScalar u = RealScalar(hi) * (RealScalar(1) / RealScalar(4294967296.0));
    return RealScalar(2) * u - RealScalar(1);
  }
};

/** \internal
 *
 * In-place LU factorization of the symmetric tridiagonal matrix \f$ T - \lambda I \f$ by
 * Gaussian elimination with partial pivoting and implicit row interchanges (LAPACK's xLAGTF):
 *   \f$ T - \lambda I = P L U \f$,
 * with \c P a permutation, \c L unit lower bidiagonal, and \c U upper triangular with up to two
 * non-zero super-diagonals. This is the careful factorization that makes inverse iteration on the
 * deliberately near-singular \f$ (T - \lambda I) \f$ (\f$ \lambda \f$ \f$ \approx \f$ a true
 * eigenvalue) robust: the pivots are never replaced here -- tiny ones are handled, when they would
 * otherwise overflow the back-solve, by tridiagonal_lagts().
 *
 * On entry \a d holds the diagonal of T, and both \a du and \a dl hold the (shared) off-diagonal of
 * the symmetric T. On exit:
 *  - \a d   holds the \c n diagonal elements of U,
 *  - \a du  holds the \c n-1 first super-diagonal elements of U,
 *  - \a du2 holds the \c n-2 second super-diagonal elements of U,
 *  - \a dl  holds the \c n-1 multipliers of L,
 *  - \a piv[k] (k < n-1) is 1 if rows k and k+1 were interchanged at step k, 0 otherwise.
 *
 * \param[in,out] d   diagonal (length \c n).
 * \param[in,out] du  first super-diagonal (length \c n-1, indices 0..n-2).
 * \param[in,out] dl  sub-diagonal / L multipliers (length \c n-1, indices 0..n-2).
 * \param[out]    du2 second super-diagonal of U (length \c n-2, indices 0..n-3).
 * \param[out]    piv row-interchange flags (length \c n).
 * \param[in]     lambda the shift \f$ \lambda \f$.
 * \param[in]     n      dimension of T.
 */
template <typename RealScalar>
void tridiagonal_lagtf(RealScalar* d, RealScalar* du, RealScalar* dl, RealScalar* du2, Index* piv, RealScalar lambda,
                       Index n) {
  d[0] -= lambda;
  piv[n - 1] = 0;
  if (n == 1) return;

  RealScalar scale1 = numext::abs(d[0]) + numext::abs(du[0]);
  for (Index k = 0; k < n - 1; ++k) {
    d[k + 1] -= lambda;
    RealScalar scale2 = numext::abs(dl[k]) + numext::abs(d[k + 1]);
    if (k < n - 2) scale2 += numext::abs(du[k + 1]);
    // Scaled magnitudes of the two candidate pivots; pick the larger (partial pivoting).
    const RealScalar piv1 = numext::is_exactly_zero(d[k]) ? RealScalar(0) : numext::abs(d[k]) / scale1;
    if (numext::is_exactly_zero(dl[k])) {
      piv[k] = 0;
      scale1 = scale2;
      if (k < n - 2) du2[k] = RealScalar(0);
    } else {
      const RealScalar piv2 = numext::abs(dl[k]) / scale2;
      if (piv2 <= piv1) {
        // No interchange: eliminate using the diagonal pivot d[k] (non-zero, since piv2 > 0).
        piv[k] = 0;
        scale1 = scale2;
        dl[k] = dl[k] / d[k];
        d[k + 1] -= dl[k] * du[k];
        if (k < n - 2) du2[k] = RealScalar(0);
      } else {
        // Interchange rows k and k+1; the sub-diagonal entry becomes the pivot.
        piv[k] = 1;
        const RealScalar mult = d[k] / dl[k];
        d[k] = dl[k];
        const RealScalar temp = d[k + 1];
        d[k + 1] = du[k] - mult * temp;
        if (k < n - 2) {
          du2[k] = du[k + 1];
          du[k + 1] = -mult * du2[k];
        }
        du[k] = temp;
        dl[k] = mult;
      }
    }
  }
}

/** \internal
 *
 * Overflow-safe solution of \f$ (T - \lambda I) x = b \f$ from the xLAGTF factorization
 * \f$ T - \lambda I = P L U \f$ (LAPACK's xLAGTS with JOB = -1). The right-hand side \a b is
 * overwritten with the solution. Diagonal elements of \c U that are so small the back-solve would
 * overflow are perturbed (by a quantity of order \f$ \mathrm{eps} \cdot \|U\| \f$, with the sign of
 * the pivot, doubling until safe) -- exactly the situation produced by inverse iteration, where
 * \f$ \lambda \f$ is intentionally a near-eigenvalue and \f$ U \f$ is nearly singular.
 *
 * \param[in]     d   diagonal of U (length \c n).
 * \param[in]     rcp elementwise reciprocal \c 1/d (length \c n), precomputed once by the caller with
 *                    an exact division and reused across the inverse-iteration solves. Entries for a
 *                    vanishing pivot may be \c inf; they are only consumed on the fast path, which is
 *                    not taken for such pivots.
 * \param[in]     du  first super-diagonal of U (length \c n-1).
 * \param[in]     dl  sub-diagonal multipliers of L (length \c n-1).
 * \param[in]     du2 second super-diagonal of U (length \c n-2).
 * \param[in]     piv row-interchange flags from tridiagonal_lagtf() (length \c n).
 * \param[in,out] b   right-hand side on input, solution on output (length \c n).
 * \param[in]     n   dimension of T.
 */
template <typename RealScalar>
void tridiagonal_lagts(const RealScalar* d, const RealScalar* rcp, const RealScalar* du, const RealScalar* dl,
                       const RealScalar* du2, const Index* piv, RealScalar* b, Index n) {
  const RealScalar eps = NumTraits<RealScalar>::epsilon();
  const RealScalar sfmin = (std::numeric_limits<RealScalar>::min)();
  const RealScalar bignum = RealScalar(1) / sfmin;

  // Perturbation floor: eps times the largest magnitude entry of U.
  using ArrayBuf = Map<const Array<RealScalar, Dynamic, 1>>;
  RealScalar tol = ArrayBuf(d, n).abs().maxCoeff();
  if (n > 1) tol = numext::maxi(tol, ArrayBuf(du, n - 1).abs().maxCoeff());
  if (n > 2) tol = numext::maxi(tol, ArrayBuf(du2, n - 2).abs().maxCoeff());
  tol *= eps;
  if (numext::is_exactly_zero(tol)) tol = eps;

  // Forward substitution: apply P then L^{-1}. Written branchlessly in the (data-dependent, roughly
  // 50/50) pivot-interchange flag so a misprediction cannot stall the serial recurrence; both arms
  // are evaluated and selected, reproducing the two branches bit for bit.
  for (Index k = 1; k < n; ++k) {
    const RealScalar dlk = dl[k - 1];
    const RealScalar bkm1 = b[k - 1];
    const RealScalar bk = b[k];
    const bool interchange = piv[k - 1] != 0;
    b[k - 1] = interchange ? bk : bkm1;
    b[k] = interchange ? (bkm1 - dlk * bk) : (bk - dlk * bkm1);
  }

  // Back substitution: solve U x = b. Away from the singularity the pivot is safe and the quotient is
  // a reciprocal multiply (rcp[k] = 1/d[k], exact and precomputed once for all the solves); only a
  // dangerously small pivot -- the rare near-singular element of the deliberately singular system --
  // falls back to the LAPACK xLAGTS perturbation/rescale with a true division. The guard is written so
  // the fast path is taken exactly when the original would divide by the unperturbed pivot, and is
  // virtually always predicted taken.
  for (Index k = n - 1; k >= 0; --k) {
    RealScalar temp = b[k];
    if (k <= n - 2) temp -= du[k] * b[k + 1];
    if (k <= n - 3) temp -= du2[k] * b[k + 2];
    const RealScalar ak = d[k];
    const RealScalar absak = numext::abs(ak);
    if (EIGEN_PREDICT_TRUE(absak >= RealScalar(1) || (absak >= sfmin && numext::abs(temp) <= absak * bignum))) {
      b[k] = temp * rcp[k];
    } else {
      // Tiny pivot: perturb (or rescale) exactly as LAPACK xLAGTS so the quotient cannot overflow.
      RealScalar a = ak;
      RealScalar t = temp;
      RealScalar pert = (a >= RealScalar(0)) ? tol : -tol;
      while (true) {
        const RealScalar aa = numext::abs(a);
        if (aa < RealScalar(1)) {
          if (aa < sfmin) {
            if (numext::is_exactly_zero(aa) || numext::abs(t) * sfmin > aa) {
              a += pert;
              pert *= RealScalar(2);
              continue;
            } else {
              t *= bignum;
              a *= bignum;
            }
          } else if (numext::abs(t) > aa * bignum) {
            a += pert;
            pert *= RealScalar(2);
            continue;
          }
        }
        break;
      }
      b[k] = t / a;
    }
  }
}

/** \internal
 *
 * Computes the eigenvectors for one contiguous block of columns [\a j_lo, \a j_hi) of a real
 * symmetric tridiagonal matrix T by inverse iteration. This is the unit of parallel work in
 * tridiagonal_inverse_iteration(): the block boundaries are aligned to cluster boundaries
 * (\c clstart[j_lo]==j_lo, and likewise \a j_hi), so the intra-cluster modified Gram-Schmidt for
 * any column in the block references only earlier columns of the same block. Two disjoint blocks
 * therefore share no state, and the result is independent of how the columns are partitioned -- in
 * particular bit-identical for any number of threads.
 *
 * The per-column shifts \a xj_scaled and cluster starts \a clstart are precomputed once by the
 * caller (see tridiagonal_inverse_iteration()) so this routine reproduces the serial sweep exactly.
 *
 * \param[in]  sdiag     normalized diagonal of T (length \c n).
 * \param[in]  ssub      normalized sub-diagonal of T (length \c n-1).
 * \param[in]  xj_scaled per-column normalized, perturbed shifts (length \c m).
 * \param[in]  clstart   per-column index of the first column of its cluster (length \c m).
 * \param[in]  n         dimension of T.
 * \param[in]  onenrm    infinity norm of the normalized T (strictly positive).
 * \param[in]  dtpcrt    inverse-iteration growth threshold.
 * \param[in]  maxits    maximum inverse-iteration steps per column.
 * \param[in]  extra     number of safety iterations past the growth threshold.
 * \param[out] eivecs    eigenvectors; only columns [\a j_lo, \a j_hi) are written (sized \c n x \c m).
 * \param[in]  j_lo,j_hi half-open range of columns to compute, aligned to cluster boundaries.
 * \returns the number of columns in the block that did not converge within \a maxits iterations.
 */
template <typename RealScalar, typename EivecType>
Index tridiagonal_inverse_iteration_block(const RealScalar* sdiag, const RealScalar* ssub, const RealScalar* xj_scaled,
                                          const Index* clstart, Index n, RealScalar onenrm, RealScalar dtpcrt,
                                          int maxits, int extra, EivecType& eivecs, Index j_lo, Index j_hi) {
  typedef Matrix<RealScalar, Dynamic, 1> RealVectorType;
  const RealScalar eps = NumTraits<RealScalar>::epsilon();

  // Work arrays for the LU factors of T - xj*I (reused across the block's columns) and the iterate.
  // lu_rcp holds the exact reciprocal of the U diagonal, computed once per factorization and reused by
  // the inverse-iteration back-solves. Being locals, each thread gets its own copy.
  RealVectorType lu_d(n), lu_dl(n), lu_du(n), lu_du2(n), lu_rcp(n), b(n);
  Matrix<Index, Dynamic, 1> piv(n);

  Index nonconv = 0;  // columns that exhausted maxits without satisfying the growth test
  for (Index j = j_lo; j < j_hi; ++j) {
    const RealScalar xj = xj_scaled[j];
    const Index gpind = clstart[j];  // first column of j's cluster (>= j_lo by block alignment)

    // Factor T - xj*I = P L U once; the inverse-iteration loop reuses it.
    lu_d = Map<const RealVectorType>(sdiag, n);
    lu_du.head(n - 1) = Map<const RealVectorType>(ssub, n - 1);
    lu_dl.head(n - 1) = Map<const RealVectorType>(ssub, n - 1);
    tridiagonal_lagtf<RealScalar>(lu_d.data(), lu_du.data(), lu_dl.data(), lu_du2.data(), piv.data(), xj, n);
    // Exact reciprocals of the U diagonal (pdiv, not the approximate preciprocal of cwiseInverse());
    // a vanishing pivot yields inf, which the back-solve fast path never consumes.
    lu_rcp.array() = RealScalar(1) / lu_d.array();

    // Deterministic pseudo-random start vector (seeded per column for thread independence).
    inverse_iteration_rng rng(numext::uint64_t(j) + 1);
    for (Index i = 0; i < n; ++i) b[i] = rng.template next<RealScalar>();

    int nrmchk = 0;
    bool converged = false;
    for (int its = 0; its < maxits; ++its) {
      // Scale the right-hand side so the near-singular solve neither overflows nor underflows.
      const RealScalar bmax = b.cwiseAbs().maxCoeff();
      if (numext::is_exactly_zero(bmax)) break;  // degenerate iterate; cannot grow -> not converged
      const RealScalar scl = RealScalar(n) * onenrm * numext::maxi(eps, numext::abs(lu_d[n - 1])) / bmax;
      b *= scl;
      tridiagonal_lagts<RealScalar>(lu_d.data(), lu_rcp.data(), lu_du.data(), lu_dl.data(), lu_du2.data(), piv.data(),
                                    b.data(), n);

      // Modified Gram-Schmidt against the already-accepted eigenvectors of this cluster.
      for (Index i = gpind; i < j; ++i) b -= b.dot(eivecs.col(i)) * eivecs.col(i);

      const RealScalar nrm = b.cwiseAbs().maxCoeff();
      if (nrm < dtpcrt) continue;          // not yet grown into the eigenvector; iterate
      if (++nrmchk < extra + 1) continue;  // a couple of safety iterations past the threshold
      converged = true;
      break;
    }
    // LAPACK xSTEIN reports vectors that never satisfied the growth test; mirror that. The column is
    // still written (best effort), but the caller surfaces the count as ComputationInfo::NoConvergence.
    if (!converged) ++nonconv;

    // Normalize to unit 2-norm with a deterministic sign (largest-magnitude entry positive). The
    // iterate can carry huge entries (~1/eps times the start) on a nearly singular solve, so divide
    // by the infinity norm first -- otherwise squaredNorm() would overflow and zero out the vector.
    Index jmax = 0;
    const RealScalar binf = b.cwiseAbs().maxCoeff(&jmax);
    if (!numext::is_exactly_zero(binf)) b /= binf;
    const RealScalar nrm2 = b.norm();
    RealScalar scl = numext::is_exactly_zero(nrm2) ? RealScalar(1) : RealScalar(1) / nrm2;
    if (b[jmax] < RealScalar(0)) scl = -scl;
    eivecs.col(j) = b * scl;
  }
  return nonconv;
}

/** \internal
 *
 * Computes eigenvectors of a real symmetric tridiagonal matrix T by inverse iteration, given a set
 * of already-computed eigenvalues (LAPACK's xSTEIN driver, built on tridiagonal_lagtf() /
 * tridiagonal_lagts()). For each requested eigenvalue \f$ \lambda_j \f$ the routine factors
 * \f$ T - \lambda_j I \f$ and applies a few steps of inverse iteration from a deterministic
 * pseudo-random start, reorthogonalizing (modified Gram-Schmidt) against the eigenvectors of any
 * tightly clustered neighbours so that a degenerate cluster yields an orthonormal basis.
 *
 * The eigenvalues are assumed sorted in non-decreasing order (as produced by spectral bisection or
 * the QR algorithm). The whole matrix is treated as a single block: an eigenvalue belonging to a
 * disconnected diagonal block (separated by a zero off-diagonal) factors to a near-singular U only
 * on its own block, so the iterate is automatically supported there; reorthogonalization across
 * blocks is harmless because those eigenvectors are already orthogonal.
 *
 * \a eivals may be an arbitrary subset of the spectrum (e.g. a bisection range): exactly one column
 * is produced per supplied eigenvalue. Reorthogonalization only ever runs among the supplied
 * eigenvalues, so the output columns are mutually orthonormal but are not orthogonalized against
 * cluster members omitted from \a eivals. A subset that splits a numerically degenerate cluster
 * therefore returns an arbitrary (still orthonormal) basis of the requested slice rather than a
 * canonical one; pass the whole cluster when that distinction matters.
 *
 * \param[in]  diag    diagonal of T (length \c n).
 * \param[in]  subdiag sub-diagonal of T (length \c n-1).
 * \param[in]  eivals  the eigenvalues whose eigenvectors are wanted, non-decreasing (length \c m).
 * \param[out] eivecs  filled with the \c m eigenvectors as its columns (must be sized \c n x \c m);
 *                     column \c j is a unit-norm eigenvector for \c eivals[j].
 * \returns the number of eigenvectors that did not converge within the inverse-iteration step limit
 *          (0 on full success); the caller maps a non-zero count to ComputationInfo::NoConvergence.
 */
template <typename DiagType, typename SubdiagType, typename EivalType, typename EivecType>
Index tridiagonal_inverse_iteration(const DiagType& diag, const SubdiagType& subdiag, const EivalType& eivals,
                                    EivecType& eivecs) {
  typedef typename DiagType::Scalar RealScalar;
  EIGEN_STATIC_ASSERT(NumTraits<RealScalar>::IsInteger == 0 && NumTraits<RealScalar>::IsComplex == 0,
                      THIS_FUNCTION_IS_NOT_FOR_INTEGER_OR_COMPLEX_TYPES)

  const Index n = diag.size();
  const Index m = eivals.size();
  if (n == 0 || m == 0) return 0;
  if (n == 1) {
    eivecs.setOnes();
    return 0;
  }

  const RealScalar eps = NumTraits<RealScalar>::epsilon();

  // Normalize T (and the shifts) to O(1) so the deliberately near-singular factor/solve cannot
  // overflow or underflow; eigenvectors are invariant under this uniform scaling. Divide each entry
  // directly by the largest magnitude rather than multiplying by its reciprocal: when that magnitude
  // is subnormal, 1/scale overflows to infinity (which would disable the normalization and let the
  // iterate underflow to an all-zero "eigenvector"), whereas entry/scale stays O(1) and finite.
  RealScalar scale = diag.cwiseAbs().maxCoeff();
  scale = numext::maxi(scale, subdiag.cwiseAbs().maxCoeff());
  if (numext::is_exactly_zero(scale)) scale = RealScalar(1);  // T == 0: any orthonormal basis works
  const Matrix<RealScalar, Dynamic, 1> sdiag = diag.array() / scale;
  const Matrix<RealScalar, Dynamic, 1> ssub = subdiag.array() / scale;

  // Infinity norm of the scaled T: max_i (|e_{i-1}| + |d_i| + |e_i|), missing boundary off-diagonals zero.
  RealScalar onenrm = numext::abs(sdiag[0]) + numext::abs(ssub[0]);
  onenrm = numext::maxi(onenrm, numext::abs(sdiag[n - 1]) + numext::abs(ssub[n - 2]));
  if (n > 2)
    onenrm = numext::maxi(onenrm, (ssub.head(n - 2).array().abs() + sdiag.segment(1, n - 2).array().abs() +
                                   ssub.segment(1, n - 2).array().abs())
                                      .maxCoeff());
  if (numext::is_exactly_zero(onenrm)) onenrm = RealScalar(1);  // T == 0: any orthonormal basis works

  // Cluster threshold and convergence threshold (LAPACK xSTEIN constants).
  const RealScalar ortol = RealScalar(1e-3) * onenrm;
  const RealScalar dtpcrt = numext::sqrt(RealScalar(0.1) / RealScalar(n));
  const int maxits = 5;
  const int extra = 2;

  // Pre-pass (sequential, O(m)): the perturbed, scaled shift xj_scaled[j] and the cluster start
  // clstart[j] (index of the first eigenvalue of j's cluster). Reproduced exactly as the serial sweep
  // would, so the parallel column blocks below match the serial result bit for bit. A shift sitting on
  // top of its predecessor is nudged up by pertol so the two factorizations stay distinct; a gap larger
  // than ortol starts a new cluster (within which the eigenvectors are reorthogonalized). The
  // perturbation is at most ~pertol << ortol, so it never moves a column across a cluster boundary.
  Matrix<RealScalar, Dynamic, 1> xj_scaled(m);
  Matrix<Index, Dynamic, 1> clstart(m);
  {
    Index gpind = 0;
    RealScalar xjm = RealScalar(0);  // previous (possibly perturbed) shift
    for (Index j = 0; j < m; ++j) {
      RealScalar xj = eivals[j] / scale;
      if (j > 0) {
        const RealScalar pertol = RealScalar(10) * numext::abs(eps * xj);
        if (xj - xjm < pertol) xj = xjm + pertol;
      }
      if (j == 0 || numext::abs(xj - xjm) > ortol) gpind = j;
      xj_scaled[j] = xj;
      clstart[j] = gpind;
      xjm = xj;
    }
  }

  const RealScalar* sdiag_p = sdiag.data();
  const RealScalar* ssub_p = ssub.data();
  const RealScalar* xj_p = xj_scaled.data();
  const Index* cl_p = clstart.data();

  // Distinct clusters are independent, so the columns split across threads with a single fork/join:
  // thread t owns a contiguous slice of columns snapped to cluster boundaries (so no cluster straddles
  // two threads) and writes only those columns of eivecs. No communication until the join, and the
  // result is bit-identical to the serial path for any thread count. Each thread tallies its own
  // non-converged columns; the reduction sums them into the value returned to the caller.
  Index nonconv = 0;
#if defined(EIGEN_HAS_OPENMP)
  int nthreads = 1;
  // Don't nest inside an existing parallel region, and only fork when there is enough work to amortize
  // the thread overhead. Each column costs ~one O(n) factorization plus a few O(n) back-solves.
  if (omp_get_num_threads() == 1) {
    // Work ~ m columns * n rows (a factorization plus a few O(n) back-solves each). kMinTaskSize is the
    // minimum such work to give each thread before adding one: m*n / kMinTaskSize threads, capped at
    // the pool size. The value is measured -- at m*n = 2*kMinTaskSize (the two-thread point) inverse
    // iteration is already ~2x faster than serial, with the gain growing to the core count for larger n.
    const double work = m * n;
    const double kMinTaskSize = 2048.0;
    const Index work_threads = Index(work / kMinTaskSize);
    nthreads = int(numext::maxi(Index(1), numext::mini(work_threads, Index(Eigen::nbThreads()))));
  }
  if (nthreads > 1) {
#pragma omp parallel num_threads(nthreads) reduction(+ : nonconv)
    {
      const Index nt = omp_get_num_threads();
      const Index tid = omp_get_thread_num();
      // Balanced split snapped up to the next cluster boundary (clstart[j]==j marks a cluster start).
      // Adjacent threads snap the same raw index to the same boundary, so the blocks tile [0, m)
      // exactly; a thread whose whole share falls inside one cluster simply gets an empty block.
      Index lo = tid * m / nt;
      Index hi = (tid + 1) * m / nt;
      while (lo < m && cl_p[lo] != lo) ++lo;
      while (hi < m && cl_p[hi] != hi) ++hi;
      nonconv += tridiagonal_inverse_iteration_block<RealScalar>(sdiag_p, ssub_p, xj_p, cl_p, n, onenrm, dtpcrt, maxits,
                                                                 extra, eivecs, lo, hi);
    }
  } else
#endif
  {
    nonconv = tridiagonal_inverse_iteration_block<RealScalar>(sdiag_p, ssub_p, xj_p, cl_p, n, onenrm, dtpcrt, maxits,
                                                              extra, eivecs, 0, m);
  }
  return nonconv;
}

/** \internal
 *
 * Rayleigh-Ritz refinement of clustered eigenvectors produced by inverse iteration. Within a
 * cluster (consecutive eigenvalues closer than \c 1e-3 * ||T||, matching the inverse-iteration
 * grouping) the computed vectors form an orthonormal basis of the invariant subspace, but they are
 * an arbitrary basis, so each \f$ (\lambda_i, v_i) \f$ pair carries a residual up to the cluster's
 * eigenvalue spread. For each such cluster of size \f$ c > 1 \f$ this forms the \f$ c \times c \f$
 * Rayleigh quotient \f$ B = V_c^\top T V_c \f$ (using the tridiagonal structure for \f$ T V_c \f$),
 * diagonalizes it, and rotates the basis into the Ritz vectors \f$ V_c Q \f$ with the (more
 * accurate) Ritz values, recovering near-machine-precision residuals.
 *
 * A cluster is refined only when its vectors are not already accurate (worst residual above a small
 * multiple of eps). A group that is well separated relative to the working precision -- the common
 * case in double precision, where \c 1e-3 * ||T|| spans many representable eigenvalues -- has
 * individual inverse-iteration vectors that are *better* than any subspace-based Ritz basis, so
 * refining it would increase the residual; such groups, singleton clusters, and a fully
 * non-degenerate spectrum are left untouched.
 *
 * The two thresholds play complementary roles, so no separate float/double tolerances are needed.
 * The grouping threshold \c 1e-3 * ||T|| is *relative* and precision-independent (the LAPACK xSTEIN
 * constant): the boundary orthogonality error is ~ \c 1e3 * eps, a fixed number of ulp in either
 * precision, and the grouping must match the driver's reorthogonalized clusters for \f$ B \f$ to be a
 * valid projection. All precision dependence lives in the eps-scaled refinement gate: that is what
 * skips a well-separated double cluster yet refines a genuinely degenerate float one.
 *
 * Only the eigenvectors are refined; the eigenvalues are left as supplied (within a degenerate
 * cluster they already agree with the Ritz values to a few ulp, and keeping them preserves the
 * contract that the returned vectors correspond to the given eigenvalues).
 *
 * \param[in]     diag    diagonal of T (length \c n).
 * \param[in]     subdiag sub-diagonal of T (length \c n-1).
 * \param[in]     eivals  the cluster eigenvalues, non-decreasing (length \c m); read, not modified.
 * \param[in,out] eivecs  the \c n x \c m eigenvectors, refined in place.
 */
template <typename DiagType, typename SubdiagType, typename EivalType, typename EivecType>
void tridiagonal_rayleigh_ritz_refine(const DiagType& diag, const SubdiagType& subdiag, const EivalType& eivals,
                                      EivecType& eivecs) {
  typedef typename DiagType::Scalar RealScalar;
  typedef Matrix<RealScalar, Dynamic, Dynamic> DenseType;
  const Index n = diag.size();
  const Index m = eivals.size();
  if (n < 2 || m < 2) return;

  // Inf-norm of T and the cluster threshold, matching tridiagonal_inverse_iteration()'s grouping.
  RealScalar onenrm = numext::abs(diag[0]) + numext::abs(subdiag[0]);
  onenrm = numext::maxi(onenrm, numext::abs(diag[n - 1]) + numext::abs(subdiag[n - 2]));
  for (Index i = 1; i < n - 1; ++i)
    onenrm = numext::maxi(onenrm, numext::abs(subdiag[i - 1]) + numext::abs(diag[i]) + numext::abs(subdiag[i]));
  const RealScalar ortol = RealScalar(1e-3) * onenrm;
  if (!(ortol > RealScalar(0))) return;  // T == 0: nothing to refine
  const RealScalar inv_onenrm = RealScalar(1) / onenrm;
  const RealScalar refine_threshold = RealScalar(16) * NumTraits<RealScalar>::epsilon();

  SelfAdjointEigenSolver<DenseType> block_solver;
  // Per-cluster scratch, hoisted out of the loop and reused across clusters: each assignment below
  // resizes only when a cluster is larger than any seen so far, so the common many-small-clusters case
  // allocates once rather than per cluster. Bsym is the distinct target for the symmetrization, which
  // also lets B + B.transpose() avoid an aliasing temporary.
  DenseType Vc, TVc, B, Bsym;
  Index s = 0;
  while (s < m) {
    Index t = s;
    while (t + 1 < m && (eivals[t + 1] - eivals[t]) < ortol) ++t;
    const Index c = t - s + 1;
    if (c > 1) {
      // Project T onto the cluster basis. T V_c via the tridiagonal structure:
      //   (T V_c)(i,:) = e_{i-1} V_c(i-1,:) + d_i V_c(i,:) + e_i V_c(i+1,:).
      Vc = eivecs.middleCols(s, c);  // copy so the write-back below cannot alias
      TVc = diag.asDiagonal() * Vc;
      TVc.topRows(n - 1).noalias() += subdiag.asDiagonal() * Vc.bottomRows(n - 1);
      TVc.bottomRows(n - 1).noalias() += subdiag.asDiagonal() * Vc.topRows(n - 1);
      // Worst per-vector residual in the cluster (scaled by 1/||T|| so it cannot overflow). Skip the
      // refinement when the cluster is already at machine precision -- refining it would only add the
      // subspace error to vectors that are individually more accurate.
      const RealScalar cluster_resid =
          ((TVc - Vc * eivals.segment(s, c).asDiagonal()) * inv_onenrm).colwise().norm().maxCoeff();
      if (cluster_resid > refine_threshold) {
        B.noalias() = Vc.transpose() * TVc;
        Bsym = B + B.transpose();  // enforce exact symmetry (distinct target avoids an aliasing temp)
        Bsym *= RealScalar(0.5);
        block_solver.compute(Bsym, ComputeEigenvectors);
        eivecs.middleCols(s, c).noalias() = Vc * block_solver.eigenvectors();
      }
    }
    s = t + 1;
  }
}

}  // namespace internal
}  // namespace Eigen

#endif  // EIGEN_TRIDIAGONAL_INVERSE_ITERATION_H
