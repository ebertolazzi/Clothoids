// Benchmark: GPU Conjugate Gradient via gpu::DeviceMatrix operators.
//
// Shows the path to running Eigen's CG on GPU with minimal code changes.
// The gpu::DeviceMatrix benchmark mirrors Eigen's conjugate_gradient() line-by-line.
// A raw cuBLAS device-pointer-mode implementation is included as a lower bound.
//
// The only change needed in Eigen's CG template to support gpu::DeviceMatrix:
//   Line 34:  typedef Dest VectorType;  (instead of Matrix<Scalar, Dynamic, 1>)
//
// Usage:
//   cmake --build build-bench-gpu --target bench_gpu_cg_sync
//   ./build-bench-gpu/bench_gpu_cg_sync
// SPDX-FileCopyrightText: The Eigen Authors
// SPDX-License-Identifier: MPL-2.0

#include <benchmark/benchmark.h>

#include <Eigen/Sparse>
#include <unsupported/Eigen/GPU>
#include <cusparse.h>

using namespace Eigen;

using Scalar = double;
using RealScalar = double;
using Vec = Matrix<Scalar, Dynamic, 1>;
using SpMat = SparseMatrix<Scalar, ColMajor, int>;

static SpMat make_spd(Index n) {
  SpMat A(n, n);
  A.reserve(VectorXi::Constant(n, 3));
  for (Index i = 0; i < n; ++i) {
    A.insert(i, i) = 4.0;
    if (i > 0) A.insert(i, i - 1) = -1.0;
    if (i < n - 1) A.insert(i, i + 1) = -1.0;
  }
  A.makeCompressed();
  return A;
}

static void cuda_warmup() {
  static bool done = false;
  if (!done) {
    void* p;
    cudaMalloc(&p, 1);
    cudaFree(p);
    done = true;
  }
}

// ==========================================================================
// GPU CG using gpu::DeviceMatrix operators — mirrors Eigen's conjugate_gradient()
// ==========================================================================
//
// Compare with Eigen/src/IterativeLinearSolvers/ConjugateGradient.h lines 29-84.
// Left column: Eigen CG code.  Right column: this benchmark.
//
//   Eigen CG                              GPU CG (this benchmark)
//   --------                              -----------------------
//   VectorType residual = rhs - mat * x;  residual.copyFrom(ctx, rhs);  [x=0 so r=b]
//   RealScalar rhsNorm2 = rhs.sqNorm();   RealScalar rhsNorm2 = rhs.squaredNorm();
//   ...
//   tmp.noalias() = mat * p;              tmp.noalias() = mat * p;  [identical]
//   Scalar alpha = absNew / p.dot(tmp);   Scalar alpha = absNew / p.dot(tmp);  [identical]
//   x += alpha * p;                       x += alpha * p;  [identical]
//   residual -= alpha * tmp;              residual -= alpha * tmp;  [identical]
//   residualNorm2 = residual.sqNorm();    residualNorm2 = residual.squaredNorm();  [identical]
//   ...
//   p = z + beta * p;                     p *= beta; p += z;  [equivalent, no alloc]

static void BM_CG_DeviceMatrixOps(benchmark::State& state) {
  cuda_warmup();
  const Index n = state.range(0);

  SpMat A = make_spd(n);
  Vec b = Vec::Random(n);

  // One shared context: SpMV + BLAS-1 on same stream, zero event overhead.
  gpu::Context ctx;
  gpu::Context::setThreadLocal(&ctx);
  gpu::SparseContext<Scalar> spmv(ctx);
  auto mat = spmv.deviceView(A);

  // Upload RHS once.
  auto rhs = gpu::DeviceMatrix<Scalar>::fromHost(b, ctx.stream());

  for (auto _ : state) {
    // --- Eigen CG lines 34-63: initialization ---
    //   typedef Dest VectorType;                               // GPU CHANGE: was Matrix<Scalar,Dynamic,1>
    //   VectorType residual = rhs - mat * x;                   // x=0, so residual = rhs
    gpu::DeviceMatrix<Scalar> x(n, 1);
    x.setZero();
    gpu::DeviceMatrix<Scalar> residual(n, 1);
    residual.copyFrom(ctx, rhs);

    //   RealScalar rhsNorm2 = rhs.squaredNorm();
    RealScalar rhsNorm2 = rhs.squaredNorm();
    if (rhsNorm2 == 0) continue;

    RealScalar tol = 1e-10;
    const RealScalar considerAsZero = (std::numeric_limits<RealScalar>::min)();
    RealScalar threshold = numext::maxi(RealScalar(tol * tol * rhsNorm2), considerAsZero);

    //   RealScalar residualNorm2 = residual.squaredNorm();
    RealScalar residualNorm2 = residual.squaredNorm();
    if (residualNorm2 < threshold) continue;

    //   VectorType p(n);
    //   p = precond.solve(residual);                           // no preconditioner: p = residual
    gpu::DeviceMatrix<Scalar> p(n, 1);
    p.copyFrom(ctx, residual);

    //   VectorType z(n), tmp(n);
    gpu::DeviceMatrix<Scalar> z(n, 1), tmp(n, 1);

    //   auto absNew = numext::real(residual.dot(p));
    //   gpu::DeviceScalar — stays on device, no sync.
    auto absNew = residual.dot(p);  // gpu::DeviceScalar, no sync

    //   while (i < maxIters) {
    Index maxIters = 200;
    Index i = 0;
    while (i < maxIters) {
      //     tmp.noalias() = mat * p;
      tmp.noalias() = mat * p;  // SpMV, device-resident

      //     auto alpha = absNew / p.dot(tmp);
      //   gpu::DeviceScalar / gpu::DeviceScalar → device kernel, no sync!
      auto alpha = absNew / p.dot(tmp);  // gpu::DeviceScalar, no sync

      //     x += alpha * p;
      //   gpu::DeviceScalar * gpu::DeviceMatrix → device-pointer axpy, no sync!
      x += alpha * p;

      //     residual -= alpha * tmp;
      residual -= alpha * tmp;  // device-pointer axpy, no sync

      //     residualNorm2 = residual.squaredNorm();
      residualNorm2 = residual.squaredNorm();  // THE one sync per iteration

      //     if (residualNorm2 < threshold) break;
      if (residualNorm2 < threshold) break;

      //     z = precond.solve(residual);
      z.copyFrom(ctx, residual);  // no preconditioner

      //     auto absOld = std::move(absNew);
      auto absOld = std::move(absNew);  // no sync, no alloc

      //     absNew = numext::real(residual.dot(z));
      absNew = residual.dot(z);  // gpu::DeviceScalar, no sync

      //     auto beta = absNew / absOld;
      //   gpu::DeviceScalar / gpu::DeviceScalar → device kernel, no sync!
      auto beta = absNew / absOld;  // gpu::DeviceScalar, no sync

      //     p = z + beta * p;
      p *= beta;  // device-pointer scal, no host sync
      p += z;

      i++;
    }
  }

  gpu::Context::setThreadLocal(nullptr);
  state.SetItemsProcessed(state.iterations() * 200);
}

BENCHMARK(BM_CG_DeviceMatrixOps)->RangeMultiplier(4)->Range(1 << 10, 1 << 20);

// ==========================================================================
// Raw cuBLAS device-pointer-mode CG (1 sync/iter) — performance lower bound
// ==========================================================================

__global__ void scalar_div_kernel(const Scalar* a, const Scalar* b, Scalar* out) { *out = *a / *b; }
__global__ void scalar_neg_kernel(const Scalar* in, Scalar* out) { *out = -(*in); }

static void BM_CG_DevicePointerMode(benchmark::State& state) {
  cuda_warmup();
  const Index n = state.range(0);
  const int maxIters = 200;

  SpMat A = make_spd(n);
  Vec b = Vec::Random(n);

  cudaStream_t stream;
  cudaStreamCreate(&stream);
  cublasHandle_t cublas;
  cublasCreate(&cublas);
  cublasSetStream(cublas, stream);

  cusparseHandle_t cusparse;
  cusparseCreate(&cusparse);
  cusparseSetStream(cusparse, stream);

  internal::DeviceBuffer d_outer((n + 1) * sizeof(int));
  internal::DeviceBuffer d_inner(A.nonZeros() * sizeof(int));
  internal::DeviceBuffer d_vals(A.nonZeros() * sizeof(Scalar));
  cudaMemcpy(d_outer.ptr, A.outerIndexPtr(), (n + 1) * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_inner.ptr, A.innerIndexPtr(), A.nonZeros() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_vals.ptr, A.valuePtr(), A.nonZeros() * sizeof(Scalar), cudaMemcpyHostToDevice);

  cusparseSpMatDescr_t matA;
  cusparseCreateCsc(&matA, n, n, A.nonZeros(), d_outer.ptr, d_inner.ptr, d_vals.ptr, CUSPARSE_INDEX_32I,
                    CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F);

  internal::DeviceBuffer d_tmp_buf(n * sizeof(Scalar));
  cusparseDnVecDescr_t tmp_x, tmp_y;
  cusparseCreateDnVec(&tmp_x, n, d_tmp_buf.ptr, CUDA_R_64F);
  cusparseCreateDnVec(&tmp_y, n, d_tmp_buf.ptr, CUDA_R_64F);
  Scalar spmv_alpha = 1.0, spmv_beta = 0.0;
  size_t ws_size = 0;
  cusparseSpMV_bufferSize(cusparse, CUSPARSE_OPERATION_NON_TRANSPOSE, &spmv_alpha, matA, tmp_x, &spmv_beta, tmp_y,
                          CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT, &ws_size);
  internal::DeviceBuffer d_workspace(ws_size);
  cusparseDestroyDnVec(tmp_x);
  cusparseDestroyDnVec(tmp_y);

  internal::DeviceBuffer d_x(n * sizeof(Scalar)), d_r(n * sizeof(Scalar));
  internal::DeviceBuffer d_p(n * sizeof(Scalar)), d_tmp(n * sizeof(Scalar));
  internal::DeviceBuffer d_b(n * sizeof(Scalar));
  internal::DeviceBuffer d_absNew(sizeof(Scalar)), d_absOld(sizeof(Scalar));
  internal::DeviceBuffer d_pdot(sizeof(Scalar)), d_alpha(sizeof(Scalar));
  internal::DeviceBuffer d_neg_alpha(sizeof(Scalar)), d_beta(sizeof(Scalar));
  internal::DeviceBuffer d_rnorm(sizeof(RealScalar));

  cudaMemcpy(d_b.ptr, b.data(), n * sizeof(Scalar), cudaMemcpyHostToDevice);

  auto spmv = [&](Scalar* x_ptr, Scalar* y_ptr) {
    cusparseDnVecDescr_t vx, vy;
    cusparseCreateDnVec(&vx, n, x_ptr, CUDA_R_64F);
    cusparseCreateDnVec(&vy, n, y_ptr, CUDA_R_64F);
    cusparseSpMV(cusparse, CUSPARSE_OPERATION_NON_TRANSPOSE, &spmv_alpha, matA, vx, &spmv_beta, vy, CUDA_R_64F,
                 CUSPARSE_SPMV_ALG_DEFAULT, d_workspace.ptr);
    cusparseDestroyDnVec(vx);
    cusparseDestroyDnVec(vy);
  };

  for (auto _ : state) {
    cudaMemsetAsync(static_cast<Scalar*>(d_x.ptr), 0, n * sizeof(Scalar), stream);
    cudaMemcpyAsync(d_r.ptr, d_b.ptr, n * sizeof(Scalar), cudaMemcpyDeviceToDevice, stream);
    cudaMemcpyAsync(d_p.ptr, d_b.ptr, n * sizeof(Scalar), cudaMemcpyDeviceToDevice, stream);

    cublasSetPointerMode(cublas, CUBLAS_POINTER_MODE_DEVICE);
    cublasDdot(cublas, n, static_cast<Scalar*>(d_r.ptr), 1, static_cast<Scalar*>(d_p.ptr), 1,
               static_cast<Scalar*>(d_absNew.ptr));

    for (int i = 0; i < maxIters; ++i) {
      spmv(static_cast<Scalar*>(d_p.ptr), static_cast<Scalar*>(d_tmp.ptr));

      cublasDdot(cublas, n, static_cast<Scalar*>(d_p.ptr), 1, static_cast<Scalar*>(d_tmp.ptr), 1,
                 static_cast<Scalar*>(d_pdot.ptr));

      scalar_div_kernel<<<1, 1, 0, stream>>>(static_cast<Scalar*>(d_absNew.ptr), static_cast<Scalar*>(d_pdot.ptr),
                                             static_cast<Scalar*>(d_alpha.ptr));
      scalar_neg_kernel<<<1, 1, 0, stream>>>(static_cast<Scalar*>(d_alpha.ptr), static_cast<Scalar*>(d_neg_alpha.ptr));

      cublasDaxpy(cublas, n, static_cast<Scalar*>(d_alpha.ptr), static_cast<Scalar*>(d_p.ptr), 1,
                  static_cast<Scalar*>(d_x.ptr), 1);
      cublasDaxpy(cublas, n, static_cast<Scalar*>(d_neg_alpha.ptr), static_cast<Scalar*>(d_tmp.ptr), 1,
                  static_cast<Scalar*>(d_r.ptr), 1);

      cublasDnrm2(cublas, n, static_cast<Scalar*>(d_r.ptr), 1, static_cast<RealScalar*>(d_rnorm.ptr));

      RealScalar rnorm;
      cudaMemcpyAsync(&rnorm, d_rnorm.ptr, sizeof(RealScalar), cudaMemcpyDeviceToHost, stream);
      cudaStreamSynchronize(stream);
      if (rnorm * rnorm < 1e-20) break;

      cudaMemcpyAsync(d_absOld.ptr, d_absNew.ptr, sizeof(Scalar), cudaMemcpyDeviceToDevice, stream);
      cublasDdot(cublas, n, static_cast<Scalar*>(d_r.ptr), 1, static_cast<Scalar*>(d_r.ptr), 1,
                 static_cast<Scalar*>(d_absNew.ptr));

      scalar_div_kernel<<<1, 1, 0, stream>>>(static_cast<Scalar*>(d_absNew.ptr), static_cast<Scalar*>(d_absOld.ptr),
                                             static_cast<Scalar*>(d_beta.ptr));

      cublasDscal(cublas, n, static_cast<Scalar*>(d_beta.ptr), static_cast<Scalar*>(d_p.ptr), 1);
      cublasSetPointerMode(cublas, CUBLAS_POINTER_MODE_HOST);
      Scalar one = 1.0;
      cublasDaxpy(cublas, n, &one, static_cast<Scalar*>(d_r.ptr), 1, static_cast<Scalar*>(d_p.ptr), 1);
      cublasSetPointerMode(cublas, CUBLAS_POINTER_MODE_DEVICE);
    }
    cudaStreamSynchronize(stream);
  }

  state.SetItemsProcessed(state.iterations() * maxIters);
  cusparseDestroySpMat(matA);
  cusparseDestroy(cusparse);
  cublasDestroy(cublas);
  cudaStreamDestroy(stream);
}

BENCHMARK(BM_CG_DevicePointerMode)->RangeMultiplier(4)->Range(1 << 10, 1 << 20);
