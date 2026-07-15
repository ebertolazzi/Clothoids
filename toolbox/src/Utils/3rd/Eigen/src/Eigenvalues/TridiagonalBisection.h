// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2026 Rasmus Munk Larsen <rmlarsen@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0

#ifndef EIGEN_TRIDIAGONAL_BISECTION_H
#define EIGEN_TRIDIAGONAL_BISECTION_H

// IWYU pragma: private
#include "./InternalHeaderCheck.h"

namespace Eigen {

/** \eigenvalues_module \ingroup Eigenvalues_Module
 *
 * \brief Selects which eigenvalues to compute in a spectral-bisection solve.
 *
 * Passed to TridiagonalEigenSolver::computeEigenvalues(). The Sturm-sequence
 * bisection algorithm can compute an arbitrary contiguous subset of the
 * spectrum, which the implicit-QR algorithm cannot. Three modes are provided:
 *  - all(): every eigenvalue (the default),
 *  - indices(il, iu): the eigenvalues with 0-based ascending indices in the
 *    half-open range [il, iu) (i.e. the (iu - il) smallest starting at il),
 *  - values(vl, vu): the eigenvalues lying in the half-open interval [vl, vu)
 *    (lower-closed, upper-open, like indices(); an eigenvalue exactly equal to
 *    vl is selected, one equal to vu is not).
 *
 * In every mode the computed eigenvalues are returned in non-decreasing order.
 */
struct EigenvalueRange {
  enum Kind { All, ByIndex, ByValue };

  Kind kind;
  Index il;
  Index iu;
  // Endpoints are held as long double so values() does not silently collapse two nearby endpoints
  // to the same number when the solver's RealScalar is long double (they are narrowed on use).
  long double vl;
  long double vu;

  /** \returns a selector for the entire spectrum. */
  static EigenvalueRange all() { return EigenvalueRange{All, 0, 0, 0.0L, 0.0L}; }

  /** \returns a selector for the eigenvalues with 0-based indices in [il, iu). */
  static EigenvalueRange indices(Index il, Index iu) { return EigenvalueRange{ByIndex, il, iu, 0.0L, 0.0L}; }

  /** \returns a selector for the eigenvalues in the half-open interval [vl, vu). */
  static EigenvalueRange values(long double vl, long double vu) { return EigenvalueRange{ByValue, 0, 0, vl, vu}; }
};

namespace internal {

/** \internal
 *
 * Evaluates Sturm sequences for a real symmetric tridiagonal matrix T at a
 * batch of shift points, using packet math to process many points in parallel.
 *
 * For a shift \c x the routine computes the number of eigenvalues of T that are
 * less than \c x via the recurrence
 *   \f$ q_0 = \alpha_0 - x,\quad q_i = (\alpha_i - \beta_{i-1}^2 / q_{i-1}) - x \f$,
 * counting how many pivots \f$ q_i \f$ are non-positive. To avoid division by
 * zero and overflow, every pivot that drops to or below \c pivmin is treated as
 * a (small) negative pivot, counted, and clamped to \c -pivmin. The same guard
 * is applied to every pivot including the first, so the scalar tail produces
 * results bit-identical to the vectorized body.
 *
 * \param[in] alpha     diagonal of T, length \c n.
 * \param[in] beta_sq   squared off-diagonal of T, \f$ \beta_i^2 \f$, length \c n-1.
 * \param[in] n         dimension of T.
 * \param[in] pivmin    smallest allowed pivot magnitude, strictly positive,
 *                      e.g. \c safemin * max_i beta_i^2 (floored to \c safemin).
 * \param[in] eval_points  the \c num_points shift values at which to evaluate.
 * \param[out] count    on output \c count[j] is the number of eigenvalues of T
 *                      less than \c eval_points[j]. Must have room for
 *                      \c num_points entries.
 */
// Evaluates the Sturm count for \c kUnroll packets of shift points at once, fully unrolled so
// the per-lane state (x, q, running count) stays in vector registers. See tridiagonal_sturm_counts().
template <int kUnroll, typename RealScalar>
EIGEN_STRONG_INLINE void tridiagonal_sturm_block(const RealScalar* alpha, const RealScalar* beta_sq, Index n,
                                                 RealScalar pivmin, const RealScalar* eval_points, RealScalar* count,
                                                 Index start) {
  typedef typename packet_traits<RealScalar>::type Packet;
  constexpr int kPacketSize = unpacket_traits<Packet>::size;
  const Packet pivmin_p = pset1<Packet>(pivmin);
  const Packet neg_pivmin_p = pset1<Packet>(-pivmin);
  const Packet one_p = pset1<Packet>(RealScalar(1));

  Packet x[kUnroll], q[kUnroll], c[kUnroll];
  const Packet alpha_0 = pset1<Packet>(alpha[0]);
  EIGEN_UNROLL_LOOP
  for (int k = 0; k < kUnroll; ++k) {
    x[k] = ploadu<Packet>(eval_points + start + k * kPacketSize);
    q[k] = psub(alpha_0, x[k]);
    const Packet mask = pcmp_le(q[k], pivmin_p);
    c[k] = pand(one_p, mask);
    q[k] = pselect(mask, pmin(q[k], neg_pivmin_p), q[k]);
  }
  for (Index i = 1; i < n; ++i) {
    const Packet alpha_i = pset1<Packet>(alpha[i]);
    const Packet beta_sq_im1 = pset1<Packet>(beta_sq[i - 1]);
    EIGEN_UNROLL_LOOP
    for (int k = 0; k < kUnroll; ++k) {
      // q = (alpha_i - beta_{i-1}^2 / q) - x.
      q[k] = psub(psub(alpha_i, pdiv(beta_sq_im1, q[k])), x[k]);
      const Packet mask = pcmp_le(q[k], pivmin_p);
      c[k] = padd(c[k], pand(one_p, mask));
      q[k] = pselect(mask, pmin(q[k], neg_pivmin_p), q[k]);
    }
  }
  EIGEN_UNROLL_LOOP
  for (int k = 0; k < kUnroll; ++k) pstoreu(count + start + k * kPacketSize, c[k]);
}

template <typename RealScalar>
void tridiagonal_sturm_counts(const RealScalar* alpha, const RealScalar* beta_sq, Index n, RealScalar pivmin,
                              const RealScalar* eval_points, RealScalar* count, Index num_points) {
  typedef typename packet_traits<RealScalar>::type Packet;
  constexpr Index kPacketSize = Index(unpacket_traits<Packet>::size);
  // num_points is a point count and so never negative; clearing the sign bit makes that explicit to
  // the compiler, which lets the division by the power-of-two kPacketSize lower to a shift (and tidies
  // the trailing cascade loops) instead of emitting signed-division sign-correction code.
  num_points &= (std::numeric_limits<Index>::max)();
  const Index full = num_points / kPacketSize;  // number of whole packets

  // Process whole packets in register-resident blocks, cascading 8 -> 4 -> 2 -> 1 so that even a
  // small batch gets enough independent recurrences to hide the pdiv latency. Eight packets is the
  // measured sweet spot on AVX2 (more would spill the per-lane state out of vector registers).
  Index p = 0;
  for (; p + 8 <= full; p += 8)
    tridiagonal_sturm_block<8>(alpha, beta_sq, n, pivmin, eval_points, count, p * kPacketSize);
  for (; p + 4 <= full; p += 4)
    tridiagonal_sturm_block<4>(alpha, beta_sq, n, pivmin, eval_points, count, p * kPacketSize);
  for (; p + 2 <= full; p += 2)
    tridiagonal_sturm_block<2>(alpha, beta_sq, n, pivmin, eval_points, count, p * kPacketSize);
  for (; p + 1 <= full; p += 1)
    tridiagonal_sturm_block<1>(alpha, beta_sq, n, pivmin, eval_points, count, p * kPacketSize);

  // Scalar tail for the remaining (< kPacketSize) points. Bit-identical to the packet path above.
  for (Index j = full * kPacketSize; j < num_points; ++j) {
    const RealScalar xj = eval_points[j];
    RealScalar qj = alpha[0] - xj;
    RealScalar cj = RealScalar(0);
    if (qj <= pivmin) {
      cj += RealScalar(1);
      qj = numext::mini(qj, -pivmin);
    }
    for (Index i = 1; i < n; ++i) {
      qj = (alpha[i] - beta_sq[i - 1] / qj) - xj;
      if (qj <= pivmin) {
        cj += RealScalar(1);
        qj = numext::mini(qj, -pivmin);
      }
    }
    count[j] = cj;
  }
}

/** \internal
 *
 * Runs the (\c t_hi - \c t_lo) independent bisections for a contiguous block of
 * target eigenvalue indices in lock-step, writing the converged shift midpoints
 * (in the internally normalized scale) to \c out[0 .. t_hi-t_lo).
 *
 * This is the unit of parallel work in tridiagonal_bisection(): the bisections of
 * two disjoint index blocks share no state, so each thread owns one block and the
 * blocks communicate only at the join. The vectorized batch evaluator keeps every
 * lane of every active bracket busy within a block.
 *
 * \param[in] alpha       normalized diagonal of T (length \c n).
 * \param[in] beta_sq     normalized squared off-diagonal (length \c n-1).
 * \param[in] n           dimension of T.
 * \param[in] pivmin      smallest allowed pivot magnitude (see tridiagonal_sturm_counts()).
 * \param[in] bracket_lo  lower end of the initial search bracket (count == 0 below it).
 * \param[in] bracket_hi  upper end of the initial search bracket (count == n above it).
 * \param[in] t_lo, t_hi  0-based target indices for this block, half-open [t_lo, t_hi).
 * \param[in] max_iters   maximum number of bisection steps.
 * \param[in] abs_tol     absolute width below which a bracket is considered converged.
 * \param[out] out        converged midpoints for indices [t_lo, t_hi); length t_hi-t_lo.
 */
template <typename RealScalar>
void tridiagonal_bisection_block(const RealScalar* alpha, const RealScalar* beta_sq, Index n, RealScalar pivmin,
                                 RealScalar bracket_lo, RealScalar bracket_hi, Index t_lo, Index t_hi, int max_iters,
                                 RealScalar abs_tol, RealScalar* out) {
  typedef Array<RealScalar, Dynamic, 1> ArrayType;
  const Index m = t_hi - t_lo;
  if (m <= 0) return;

  // The eigenvalue with 0-based index i is the value where count(x) crosses from <= i to > i.
  ArrayType lower = ArrayType::Constant(m, bracket_lo);
  ArrayType upper = ArrayType::Constant(m, bracket_hi);
  const ArrayType targets = ArrayType::LinSpaced(m, RealScalar(t_lo), RealScalar(t_hi - 1));
  ArrayType counts(m);
  ArrayType mid = RealScalar(0.5) * (lower + upper);

  // Each eigenvalue is recorded the iteration it first converges, using a per-element criterion that
  // depends only on that element's own bracket. This makes the result independent of how the index
  // range is grouped -- in particular bitwise identical for any number of threads -- because element
  // i's converged value never depends on when its neighbors finish. (A single global convergence test
  // would instead keep refining an already-converged eigenvalue until the slowest one in its group
  // caught up, so the last bits would shift with the thread count.)
  ArrayType result = mid;
  Array<bool, Dynamic, 1> done = Array<bool, Dynamic, 1>::Constant(m, false);

  // In the early iterations many eigenvalues still share a bracket, so their (sorted) midpoints are
  // identical and the Sturm count need only be evaluated at the distinct values and scattered back.
  // Once every midpoint is distinct the brackets only ever get finer, so we stop deduplicating and
  // evaluate all midpoints directly, avoiding the per-iteration bookkeeping (and its unpredictable
  // branch). Deduplication never changes the result: count(mid[i]) is unchanged. It only pays off
  // when there are many more points than packets (otherwise the per-point kernel cost is too small
  // to offset the bookkeeping), so the scratch is allocated only then.
  constexpr int kPacketSize = unpacket_traits<typename packet_traits<RealScalar>::type>::size;
  bool deduplicate = (m >= 64 * Index(kPacketSize));
  ArrayType distinct, counts_distinct;
  if (deduplicate) {
    distinct.resize(m);
    counts_distinct.resize(m);
  }
  for (int iter = 0; iter < max_iters; ++iter) {
    if (deduplicate) {
      const RealScalar* midp = mid.data();
      RealScalar* distp = distinct.data();
      Index nd = 0;
      distp[nd++] = midp[0];
      for (Index i = 1; i < m; ++i)
        if (midp[i] != midp[i - 1]) distp[nd++] = midp[i];
      tridiagonal_sturm_counts<RealScalar>(alpha, beta_sq, n, pivmin, distp, counts_distinct.data(), nd);
      const RealScalar* cdp = counts_distinct.data();
      RealScalar* cp = counts.data();
      Index g = 0;
      cp[0] = cdp[g];
      for (Index i = 1; i < m; ++i) {
        if (midp[i] != midp[i - 1]) ++g;
        cp[i] = cdp[g];
      }
      if (nd == m) deduplicate = false;  // all distinct: finer brackets stay distinct
    } else {
      tridiagonal_sturm_counts<RealScalar>(alpha, beta_sq, n, pivmin, mid.data(), counts.data(), m);
    }
    // count(mid) <= target  =>  eigenvalue is >= mid, raise the lower bound;
    // otherwise the eigenvalue is < mid, so lower the upper bound.
    const auto raise = (counts <= targets);
    lower = raise.select(mid, lower);
    upper = raise.select(upper, mid);
    const ArrayType new_mid = RealScalar(0.5) * (lower + upper);
    // Freeze each eigenvalue at the midpoint of its bracket the first time that bracket is tight
    // enough (width within abs_tol) or stops moving. Already-frozen entries keep their value.
    const auto converged = (new_mid == mid) || ((upper - lower) <= abs_tol);
    result = (!done && converged).select(new_mid, result);
    done = done || converged;
    mid = new_mid;
    if (done.all()) break;
  }
  // Any eigenvalue that did not converge within max_iters keeps its last midpoint, written straight
  // into the caller's output (no aliasing temporary: select is coefficient-wise and out is disjoint).
  Eigen::Map<ArrayType>(out, m) = done.select(result, mid);
}

/** \internal
 *
 * Returns the number of eigenvalues of the (normalized) symmetric tridiagonal matrix T that are
 * strictly less than \c x, evaluated by the same pivot recurrence as tridiagonal_sturm_counts() but
 * counting only strictly-negative pivots. This translates the endpoints of a half-open value range
 * [vl, vu) into eigenvalue indices [count(vl), count(vu)): counting strictly-below makes an
 * eigenvalue that lands exactly on an endpoint count as not-below, so it is kept at the closed lower
 * end and dropped at the open upper end. (The batch evaluator above instead counts a pivot at the
 * floor as below, which on its own would yield (vl, vu]; the two differ only when an endpoint
 * coincides with an eigenvalue to within the pivot floor \c pivmin, where the choice is anyway
 * numerically ambiguous.)
 *
 * As in the batch evaluator, a pivot whose magnitude reaches the floor \c pivmin is replaced by
 * +/- pivmin, keeping its sign, so the following division cannot overflow.
 *
 * \param[in] alpha    normalized diagonal of T (length \c n).
 * \param[in] beta_sq  normalized squared off-diagonal (length \c n-1).
 * \param[in] n        dimension of T.
 * \param[in] pivmin   smallest allowed pivot magnitude (see tridiagonal_sturm_counts()).
 * \param[in] x        the shift at which to count.
 */
template <typename RealScalar>
Index tridiagonal_sturm_count_below(const RealScalar* alpha, const RealScalar* beta_sq, Index n, RealScalar pivmin,
                                    RealScalar x) {
  RealScalar q = alpha[0] - x;
  Index count = (q < RealScalar(0)) ? 1 : 0;
  if (numext::abs(q) < pivmin) q = numext::copysign(pivmin, q);
  for (Index i = 1; i < n; ++i) {
    q = (alpha[i] - beta_sq[i - 1] / q) - x;
    if (q < RealScalar(0)) ++count;
    if (numext::abs(q) < pivmin) q = numext::copysign(pivmin, q);
  }
  return count;
}

/** \internal
 *
 * Computes a subset of the eigenvalues of a real symmetric tridiagonal matrix T
 * by Sturm-sequence spectral bisection (cf. LAPACK's xSTEBZ), accelerated by the
 * vectorized batch evaluator tridiagonal_sturm_counts().
 *
 * The eigenvalues are bracketed using the Gershgorin disc theorem, then refined
 * by running an independent bisection per requested eigenvalue, with all active
 * brackets evaluated together in each step so the packet evaluator stays full.
 *
 * \param[in] diag     diagonal of T (length \c n).
 * \param[in] subdiag  sub-diagonal of T (length \c n-1).
 * \param[in] range    which eigenvalues to compute (see EigenvalueRange).
 * \param[in] abs_tol  requested absolute accuracy; the effective tolerance is
 *                     max(abs_tol, eps * ||T||). Pass 0 for full precision.
 * \param[out] eivalues  filled with the selected eigenvalues in non-decreasing
 *                     order; resized to the number of selected eigenvalues.
 * \returns the number of eigenvalues written to \c eivalues.
 */
template <typename DiagType, typename SubdiagType, typename EivalType>
Index tridiagonal_bisection(const DiagType& diag, const SubdiagType& subdiag, const EigenvalueRange& range,
                            typename DiagType::Scalar abs_tol, EivalType& eivalues) {
  typedef typename DiagType::Scalar RealScalar;
  typedef Array<RealScalar, Dynamic, 1> ArrayType;
  EIGEN_STATIC_ASSERT(NumTraits<RealScalar>::IsInteger == 0 && NumTraits<RealScalar>::IsComplex == 0,
                      THIS_FUNCTION_IS_NOT_FOR_INTEGER_OR_COMPLEX_TYPES)

  const Index n = diag.size();
  if (n == 0) {
    eivalues.derived().resize(0);
    return 0;
  }

  // Normalize the matrix to O(1) to avoid overflow/underflow when squaring the
  // off-diagonal and during the Sturm recurrence; eigenvalues scale linearly, so the
  // scaling is undone at the very end. (The caller has already verified the input is
  // finite.) This mirrors the uniform scaling done in SelfAdjointEigenSolver::compute().
  // Divide each entry directly by the scale rather than multiplying by its reciprocal: when the scale
  // is subnormal, 1/scale overflows to infinity, the normalization silently disables itself, and the
  // Sturm recurrence underflows and returns wrong eigenvalues.
  RealScalar scale = diag.cwiseAbs().maxCoeff();
  if (n >= 2) scale = numext::maxi(scale, subdiag.cwiseAbs().maxCoeff());
  if (numext::is_exactly_zero(scale)) scale = RealScalar(1);

  // Local contiguous copies of the scaled matrix data, |off-diagonal|, and its square.
  const ArrayType alpha = diag.array() / scale;
  const ArrayType beta_abs = (n >= 2) ? ArrayType(subdiag.array().abs() / scale) : ArrayType(0);
  const ArrayType beta_sq = (n >= 2) ? ArrayType(beta_abs.square()) : ArrayType(0);

  // Smallest pivot allowed during the Sturm recurrence (positive, floored so
  // that a matrix with an all-zero off-diagonal still has pivmin > 0).
  const RealScalar eps = NumTraits<RealScalar>::epsilon();
  const RealScalar safemin = numext::maxi(RealScalar(1) / NumTraits<RealScalar>::highest(),
                                          (RealScalar(1) + eps) * (std::numeric_limits<RealScalar>::min)());
  const RealScalar max_beta_sq = (n >= 2) ? beta_sq.maxCoeff() : RealScalar(0);
  const RealScalar pivmin = safemin * numext::maxi(max_beta_sq, RealScalar(1));

  // Gershgorin bounds: row k has radius |beta_{k-1}| + |beta_k| (with the
  // missing boundary off-diagonals taken as zero).
  ArrayType radius(n);
  if (n == 1) {
    radius(0) = RealScalar(0);
  } else {
    radius(0) = beta_abs(0);
    radius(n - 1) = beta_abs(n - 2);
    if (n > 2) radius.segment(1, n - 2) = beta_abs.head(n - 2) + beta_abs.segment(1, n - 2);
  }
  RealScalar lambda_min = (alpha - radius).minCoeff();
  RealScalar lambda_max = (alpha + radius).maxCoeff();

  // Effective tolerance and outward expansion of the bracket so that
  // count(lambda_min) == 0 and count(lambda_max) == n (cf. LAPACK xSTEBZ).
  const RealScalar tnorm = numext::maxi(numext::abs(lambda_min), numext::abs(lambda_max));
  // Convergence tolerance: refine each bracket to ~1 ulp of ||T||. This is the role of xSTEBZ's
  // relative tolerance (it stops once b - a < RELFAC * ulp * max(|a|, |b|), RELFAC = 2), here
  // specialized to the matrix norm tnorm and folded together with any absolute tolerance the caller
  // requested.
  abs_tol = numext::maxi(abs_tol, eps * tnorm);
  // Widen the Gershgorin bracket so that, despite rounding in the Sturm recurrence, count() really does
  // reach 0 at lambda_min and n at lambda_max. The n*eps*tnorm term bounds the worst-case count error
  // accumulated over the n recurrence steps; the 2*pivmin term covers the pivot floor. The 2.1 prefactor
  // is xSTEBZ's FUDGE factor: ideally 1 would suffice, but it is taken slightly larger to stay robust on
  // sloppy arithmetic, and (per xSTEBZ) widening the bracket this way only loosens the initial search
  // interval -- it has no effect on the accuracy of the converged eigenvalues.
  const RealScalar expand = RealScalar(2.1) * (RealScalar(n) * eps * tnorm + RealScalar(2) * pivmin);
  lambda_min -= expand;
  lambda_max += expand;

  // Determine the target indices [t_lo, t_hi) and the search bracket.
  Index t_lo = 0, t_hi = n;
  RealScalar bracket_lo = lambda_min, bracket_hi = lambda_max;
  if (range.kind == EigenvalueRange::ByIndex) {
    eigen_assert(range.il >= 0 && range.il <= range.iu && range.iu <= n && "invalid eigenvalue index range");
    t_lo = range.il;
    t_hi = range.iu;
  } else if (range.kind == EigenvalueRange::ByValue) {
    eigen_assert(range.vl <= range.vu && "invalid eigenvalue value range");
    const RealScalar vl = RealScalar(range.vl) / scale;
    const RealScalar vu = RealScalar(range.vu) / scale;
    // The eigenvalues in [vl, vu) are exactly those with indices in [t_lo, t_hi), where t_lo and t_hi
    // count the eigenvalues strictly below each endpoint. Counting strictly-below puts an eigenvalue
    // that coincides with an endpoint inside the interval at the closed lower end and outside it at the
    // open upper end (see tridiagonal_sturm_count_below()).
    t_lo = tridiagonal_sturm_count_below<RealScalar>(alpha.data(), beta_sq.data(), n, pivmin, vl);
    t_hi = tridiagonal_sturm_count_below<RealScalar>(alpha.data(), beta_sq.data(), n, pivmin, vu);
    bracket_lo = numext::maxi(lambda_min, vl);
    bracket_hi = numext::mini(lambda_max, vu);
  }

  const Index m = t_hi - t_lo;
  eivalues.derived().resize(m);
  if (m <= 0) return 0;

  const int max_iters = NumTraits<RealScalar>::digits() + 2;
  ArrayType mid_all(m);
  RealScalar* out = mid_all.data();
  const RealScalar* alpha_p = alpha.data();
  const RealScalar* beta_sq_p = beta_sq.data();

  // The m eigenvalues are bisected independently, so the spectrum splits cleanly across threads with
  // a single fork/join: thread t owns the contiguous index block [t_lo + lo, t_lo + hi) and writes
  // its converged midpoints into the disjoint slice out[lo, hi). No communication until the join.
#if defined(EIGEN_HAS_OPENMP)
  int nthreads = 1;
  // Don't nest inside an existing parallel region, and only fork when there is enough work to
  // amortize the thread overhead while still leaving each thread a SIMD-friendly chunk of points.
  if (omp_get_num_threads() == 1) {
    constexpr Index kPacketSize = Index(unpacket_traits<typename packet_traits<RealScalar>::type>::size);
    // One work unit ~ one Sturm step (a packet division); kMinTaskSize is the minimum per thread.
    const double work = m * n * max_iters;
    const double kMinTaskSize = 131072.0;
    const Index work_threads = Index(work / kMinTaskSize);
    const Index point_threads = m / (8 * kPacketSize);
    const Index pb = numext::maxi(Index(1), numext::mini(work_threads, point_threads));
    nthreads = int(numext::mini(pb, Index(Eigen::nbThreads())));
  }
  if (nthreads > 1) {
#pragma omp parallel num_threads(nthreads)
    {
      const Index nt = omp_get_num_threads();
      const Index tid = omp_get_thread_num();
      // Balanced split: every thread gets floor(m/nt) or ceil(m/nt) consecutive indices, none empty.
      const Index lo = tid * m / nt;
      const Index hi = (tid + 1) * m / nt;
      tridiagonal_bisection_block<RealScalar>(alpha_p, beta_sq_p, n, pivmin, bracket_lo, bracket_hi, t_lo + lo,
                                              t_lo + hi, max_iters, abs_tol, out + lo);
    }
  } else
#endif
  {
    tridiagonal_bisection_block<RealScalar>(alpha_p, beta_sq_p, n, pivmin, bracket_lo, bracket_hi, t_lo, t_hi,
                                            max_iters, abs_tol, out);
  }

  // Undo the normalization.
  eivalues = (mid_all * scale).matrix();
  return m;
}

}  // namespace internal
}  // namespace Eigen

#endif  // EIGEN_TRIDIAGONAL_BISECTION_H
