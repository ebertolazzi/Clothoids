# Eigen GPU Module (`unsupported/Eigen/GPU`)

GPU-accelerated linear algebra for Eigen users, dispatching to NVIDIA CUDA
Math Libraries (cuBLAS, cuSOLVER, cuFFT, cuSPARSE, cuDSS). Requires CUDA 11.4+;
cuDSS features require CUDA 12.0+ and a separate cuDSS install. Header-only.

## Why this module

Eigen is the linear algebra foundation for a large ecosystem of C++ projects
in robotics (ROS, Drake, MoveIt, Pinocchio), computer vision (OpenCV, COLMAP,
Open3D), scientific computing (Ceres, Stan), and beyond. Many of these
projects run on GPU-equipped hardware but cannot use GPUs for Eigen operations
without dropping down to raw CUDA library APIs.

GPU sparse solvers are a particularly acute gap. Sparse factorization is the
bottleneck in SLAM, bundle adjustment, FEM, and nonlinear optimization --
exactly the workloads where GPU acceleration matters most. Downstream projects
like [Ceres](https://github.com/ceres-solver/ceres-solver/issues/1151) and
[COLMAP](https://github.com/colmap/colmap/issues/4018) have open requests for
GPU-accelerated sparse solvers, and third-party projects like
[cholespy](https://github.com/rgl-epfl/cholespy) exist specifically because
Eigen lacks them. The `unsupported/Eigen/GPU` module provides GPU sparse Cholesky, LDL^T,
and LU factorization via cuDSS, alongside dense solvers (cuSOLVER), matrix
products (cuBLAS), FFT (cuFFT), and sparse matrix-vector products (cuSPARSE).

Existing Eigen users should be able to move performance-critical dense or
sparse linear algebra to the GPU with minimal code changes and without
learning CUDA library APIs directly.

## Design philosophy

**CPU and GPU coexist.** There is no global compile-time switch that replaces
CPU implementations (unlike `EIGEN_USE_LAPACKE`). Users choose GPU solvers
explicitly -- `gpu::LLT<double>` vs `Eigen::LLT<MatrixXd>`,
`gpu::SparseLLT<double>` vs `SimplicialLLT<SparseMatrix<double>>` -- and both
coexist in the same binary.
This also lets users keep the factored matrix on device across multiple solves,
something impossible with compile-time replacement.

**Familiar syntax.** GPU operations use the same expression patterns as CPU
Eigen. Here is a side-by-side comparison:

```cpp
// ---- CPU (Eigen) ----               // ---- GPU (unsupported/Eigen/GPU) ----
#include <Eigen/Dense>                  #define EIGEN_USE_GPU
                                        #include <unsupported/Eigen/GPU>

// Dense
MatrixXd A = ...;                       auto d_A = gpu::DeviceMatrix<double>::fromHost(A);
MatrixXd B = ...;                       auto d_B = gpu::DeviceMatrix<double>::fromHost(B);

MatrixXd C = A * B;                     gpu::DeviceMatrix<double> d_C = d_A * d_B;
MatrixXd X = A.llt().solve(B);          gpu::DeviceMatrix<double> d_X = d_A.llt().solve(d_B);

                                        MatrixXd X = d_X.toHost();

// Sparse (using SpMat = SparseMatrix<double>)
SimplicialLLT<SpMat> llt(A);            gpu::SparseLLT<double> llt(A);
VectorXd x = llt.solve(b);              VectorXd x = llt.solve(b);
```

The GPU version reads like CPU Eigen with explicit upload/download for dense
operations, and an almost identical API for sparse solvers. Unsupported
expressions are compile errors.

**Standalone module.** `unsupported/Eigen/GPU` does not modify or depend on Eigen's Core
expression template system (`MatrixBase`, `CwiseBinaryOp`, etc.).
`DeviceMatrix` is not an Eigen expression type and does not inherit from
`MatrixBase`. The expression layer is a thin compile-time dispatch where every
supported expression maps to a single NVIDIA library call. There is no
coefficient-level evaluation, lazy fusion, or packet operations.

**Interoperability where useful.** `DeviceMatrix` provides the same operator
signatures as `Matrix` for common vector operations: `+=`, `-=`, `*=`,
`dot()`, `squaredNorm()`, `norm()`, `setZero()`, and `noalias()`. This makes
`DeviceMatrix` usable as a drop-in `VectorType` in Eigen algorithm templates
that rely on these operations. For example, Eigen's `conjugate_gradient()`
template works with `DeviceMatrix` with a single typedef change -- no
modifications to the algorithm or the expression template system. Conjugate
gradient is just the motivating example; we are open to expanding operator
coverage as needed to support other high-level Eigen algorithms on the GPU.

**Explicit over implicit.** Host-device transfers, stream management, and
library handle lifetimes are visible in the API. There are no hidden
allocations or synchronizations except where documented (e.g., `toHost()` must
synchronize to deliver data to the host).

## Key concepts

### `gpu::DeviceMatrix<Scalar>`

A typed RAII wrapper for a dense column-major matrix in GPU device memory.
This is the GPU counterpart of Eigen's `MatrixX<Scalar>`. A vector is simply
a `DeviceMatrix` with one column. All public GPU classes live in `namespace
Eigen::gpu`.

```cpp
// Upload from host
auto d_A = gpu::DeviceMatrix<double>::fromHost(A);

// Allocate uninitialized
gpu::DeviceMatrix<double> d_C(m, n);

// Download to host
MatrixXd C = d_C.toHost();

// Async download (returns a future)
auto transfer = d_C.toHostAsync();
// ... do other work ...
MatrixXd C = transfer.get();
```

`DeviceMatrix` supports expression methods that mirror Eigen's API:
`adjoint()`, `transpose()`, `triangularView<UpLo>()`,
`selfadjointView<UpLo>()`, `llt()`, `lu()`. These return lightweight
expression objects that are evaluated when assigned.

For BLAS Level-1 operations, `DeviceMatrix` also provides `dot()`, `norm()`,
`squaredNorm()`, `setZero()`, `noalias()`, and arithmetic operators
(`+=`, `-=`, `*=`) that dispatch to cuBLAS `axpy`, `nrm2`, `dot`, `scal`,
and `geam`. These are the operations needed by iterative solvers.

### `gpu::DeviceScalar<Scalar>`

A device-resident scalar value. Reductions like `dot()`, `norm()`, and
`squaredNorm()` return `DeviceScalar` instead of a host scalar, deferring
the host synchronization until the value is actually needed:

```cpp
auto dot_val = d_x.dot(d_y);          // DeviceScalar -- no sync
auto norm_sq = d_r.squaredNorm();      // DeviceScalar -- no sync
Scalar alpha = dot_val / norm_sq;      // sync here (implicit conversion)
d_x += alpha * d_p;                    // host scalar * DeviceMatrix (axpy)
```

Division between `DeviceScalar` values (real types only) is performed on
device via NPP, avoiding extra synchronizations. Small device allocations
(including `DeviceScalar`) are recycled through a thread-local
`DeviceBufferPool` to avoid `cudaMalloc`/`cudaFree` overhead in tight loops.
Pool contract: recycled pointers are safe for same-thread, same-stream reuse
(the typical iterative-solver pattern, where one `gpu::Context` drives all
work). Mixing pooled buffers across threads or CUDA streams without external
synchronization is not supported.

### `gpu::Context`

Every GPU operation needs a CUDA stream and library handles (cuBLAS eagerly,
cuSOLVER / cuBLASLt / cuSPARSE lazily on first use). `gpu::Context` bundles
these together. A single `Context` is not thread-safe -- use one per thread
(or external synchronization), since the underlying NVIDIA library handles
are not thread-safe per handle.

For simple usage, you don't need to create one -- a per-thread default context
is created lazily on first use:

```cpp
// These use the thread-local default context automatically
d_C = d_A * d_B;
d_X = d_A.llt().solve(d_B);
```

For concurrent multi-stream execution, create explicit contexts:

```cpp
gpu::Context ctx1, ctx2;
d_C1.device(ctx1) = d_A1 * d_B1;   // runs on stream 1
d_C2.device(ctx2) = d_A2 * d_B2;   // runs on stream 2 (concurrently)
```

To integrate with existing CUDA code, borrow an existing stream:

```cpp
gpu::Context ctx(my_existing_stream);  // wraps stream, does not take ownership
```

To override the thread-local default (e.g., in CG where all ops share one
context):

```cpp
gpu::Context ctx;
gpu::Context::setThreadLocal(&ctx);    // all threadLocal() calls return ctx
// ... GPU operations ...
gpu::Context::setThreadLocal(nullptr); // restore lazy-created default
```

### Linking {#eigen_gpu_linking}

The module is header-only, but each feature pulls in the corresponding NVIDIA
library at link time. cuSOLVER, cuBLASLt, and cuSPARSE are created lazily on
first use, so a translation unit that only uses cuBLAS or cuFFT does not need
to link the others:

| Feature                                 | Link flags                |
|-----------------------------------------|---------------------------|
| `DeviceMatrix`, GEMM, TRSM, SYMM, SYRK  | `-lcublas -lcublasLt`     |
| Dense solvers (LLT, LU, QR, SVD, EVD)   | `-lcusolver -lcublas`     |
| FFT (`gpu::FFT`)                        | `-lcufft -lcublas`        |
| SpMV / SpMM (`gpu::SparseContext`)      | `-lcusparse -lcublas`     |
| Sparse direct solvers (cuDSS)           | `-lcudss -lcublas`        |

cuBLAS is required by `DeviceMatrix` itself (every `Context` creates a cuBLAS
handle eagerly) and is also a runtime dependency of cuDSS, so it is the one
constant. cuDSS additionally requires `EIGEN_CUDSS` to be defined before
including `unsupported/Eigen/GPU`.

## Usage

### Matrix operations (cuBLAS)

```cpp
auto d_A = gpu::DeviceMatrix<double>::fromHost(A);
auto d_B = gpu::DeviceMatrix<double>::fromHost(B);

// GEMM: C = A * B, C = A^H * B, C = A * B^T, ...
gpu::DeviceMatrix<double> d_C = d_A * d_B;
d_C = d_A.adjoint() * d_B;
d_C = d_A * d_B.transpose();

// Scaled and accumulated
d_C += 2.0 * d_A * d_B;             // alpha=2, beta=1
d_C.device(ctx) -= d_A * d_B;       // alpha=-1, beta=1 (GEMM requires explicit context for -=)

// Triangular solve (TRSM)
d_X = d_A.triangularView<Lower>().solve(d_B);

// Symmetric/Hermitian multiply (SYMM/HEMM)
d_C = d_A.selfadjointView<Lower>() * d_B;

// Rank-k update (SYRK/HERK)
d_C.selfadjointView<Lower>().rankUpdate(d_A);  // C += A * A^H
```

### BLAS Level-1 operations

```cpp
// Dot product and norms (return DeviceScalar -- no sync until read)
auto dot_val = d_x.dot(d_y);          // cublasDdot / cublasCdotc
auto norm_val = d_r.norm();            // cublasDnrm2
double n = norm_val;                   // implicit conversion triggers sync

// Vector arithmetic (cuBLAS axpy / geam)
d_x += alpha * d_p;                    // axpy: x = x + alpha * p
d_x -= alpha * d_p;                    // axpy: x = x - alpha * p
d_x *= alpha;                          // scal: x = alpha * x
d_r.setZero();                         // cudaMemsetAsync

// DeviceScalar arithmetic (stays on device, real types only)
auto alpha = absNew / dot_val;         // device-side division via NPP
d_x += alpha * d_p;                    // DeviceScalar * DeviceMatrix (axpy with device pointer)

// Matrix add/subtract (cuBLAS geam)
gpu::DeviceMatrix<double> d_C = d_A + d_B;       // C = A + B
d_C = d_A + 2.0 * d_B;                          // C = A + 2*B
d_C = d_A - d_B;                                 // C = A - B
```

### Dense solvers (cuSOLVER)

**One-shot expression syntax** -- Convenient, re-factorizes each time:

```cpp
// Cholesky solve (potrf + potrs)
d_X = d_A.llt().solve(d_B);

// LU solve (getrf + getrs)
d_Y = d_A.lu().solve(d_B);
```

**Cached factorization** -- Factor once, solve many times:

```cpp
gpu::LLT<double> llt;
llt.compute(d_A);                    // factorize (async)
if (llt.info() != Success) { ... }   // lazy sync on first info() call
auto d_X1 = llt.solve(d_B1);        // reuses factor (async)
auto d_X2 = llt.solve(d_B2);        // reuses factor (async)
MatrixXd X2 = d_X2.toHost();

// LU with transpose solve
gpu::LU<double> lu;
lu.compute(d_A);
auto d_Y = lu.solve(d_B, gpu::GpuOp::Trans);           // A^T Y = B

// QR solve (overdetermined least squares)
gpu::QR<double> qr;
qr.compute(d_A);                     // factorize on device (async)
auto d_X = qr.solve(d_B);           // Q^H * B via ormqr, then trsm on R
MatrixXd X = d_X.toHost();

// SVD (results downloaded on access)
gpu::SVD<double> svd;
svd.compute(d_A, ComputeThinU | ComputeThinV);
VectorXd S = svd.singularValues();   // downloads to host
MatrixXd U = svd.matrixU();          // downloads to host
MatrixXd V = svd.matrixV();          // V (matches JacobiSVD)
MatrixXd VT = svd.matrixVT();        // V^T (matches cuSOLVER)

// SVD: device-side views (no D2H transfer; svd must outlive the views)
auto d_S = svd.d_singularValues();   // DeviceMatrix view of singular values
auto d_U = svd.d_matrixU();          // DeviceMatrix view of U
auto d_VT = svd.d_matrixVT();        // DeviceMatrix view of V^T

// Self-adjoint eigenvalue decomposition
gpu::SelfAdjointEigenSolver<double> es;
es.compute(d_A);
VectorXd eigenvals = es.eigenvalues();    // downloads to host
MatrixXd eigenvecs = es.eigenvectors();   // downloads to host
auto d_W = es.d_eigenvalues();            // DeviceMatrix view of eigenvalues
auto d_V = es.d_eigenvectors();           // DeviceMatrix view of eigenvectors
```

The cached API keeps the factored matrix on device, avoiding redundant
host-device transfers and re-factorizations. All solvers also accept host
matrices directly as a convenience (e.g., `gpu::LLT<double> llt(A)` or
`qr.solve(B)`), which handles upload/download internally. The `d_*` accessors
on `gpu::SVD` and `gpu::SelfAdjointEigenSolver` return non-owning
`DeviceMatrix` views so downstream cuBLAS/cuSOLVER work can chain without
round-tripping through host memory.

### Sparse direct solvers (cuDSS)

Requires cuDSS (separate install, CUDA 12.0+). Define `EIGEN_CUDSS` before
including `unsupported/Eigen/GPU`; see [Linking](#eigen_gpu_linking) for link flags.

```cpp
SparseMatrix<double> A = ...;  // symmetric positive definite
VectorXd b = ...;

// Sparse Cholesky -- one-liner
gpu::SparseLLT<double> llt(A);
VectorXd x = llt.solve(b);

// Three-phase workflow for repeated solves with the same sparsity pattern
gpu::SparseLLT<double> llt;
llt.analyzePattern(A);               // symbolic analysis (once)
llt.factorize(A);                    // numeric factorization
VectorXd x = llt.solve(b);
llt.factorize(A_new_values);         // refactorize (reuses symbolic analysis)
VectorXd x2 = llt.solve(b);

// Sparse LDL^T (symmetric indefinite)
gpu::SparseLDLT<double> ldlt(A);
VectorXd x = ldlt.solve(b);

// Sparse LU (general non-symmetric)
gpu::SparseLU<double> lu(A);
VectorXd x = lu.solve(b);
```

### FFT (cuFFT)

```cpp
gpu::FFT<float> fft;                // shares stream + cuBLAS with the
                                    // thread-local default Context
gpu::Context ctx;
gpu::FFT<float> fft_on_ctx(ctx);    // share stream + cuBLAS with an
                                    // explicit Context (e.g. for
                                    // multi-stream pipelines)

// 1D complex-to-complex
VectorXcf X = fft.fwd(x);           // forward
VectorXcf y = fft.inv(X);           // inverse (scaled by 1/n)

// 1D real-to-complex / complex-to-real
VectorXcf R = fft.fwd(r);           // returns n/2+1 complex (half-spectrum)
VectorXf  s = fft.invReal(R, n);    // C2R inverse, caller specifies n

// 2D complex-to-complex
MatrixXcf B = fft.fwd2(A);         // 2D forward
MatrixXcf C = fft.inv2(B);         // 2D inverse (scaled by 1/(rows*cols))

// Plans are cached and reused across calls with the same size/type.
```

### Sparse matrix-vector multiply (cuSPARSE)

```cpp
SparseMatrix<double> A = ...;
VectorXd x = ...;

// Host vectors (upload/download handled internally)
gpu::SparseContext<double> spmv;
VectorXd y = spmv.multiply(A, x);           // y = A * x
VectorXd z = spmv.multiplyT(A, x);          // z = A^T * x
spmv.multiply(A, x, y, 2.0, 1.0);           // y = 2*A*x + y
spmv.multiply(A, x, y, 1.0, 0.0,            // y = A^H * x (Hermitian SpMV)
              gpu::GpuOp::ConjTrans);

// Multiple RHS (SpMM)
MatrixXd Y = spmv.multiplyMat(A, X);                      // Y = A * X
MatrixXd Z = spmv.multiplyMat(A, X, gpu::GpuOp::Trans);   // Z = A^T * X

// Device-resident SpMV (sparse matrix cached on device)
gpu::Context ctx;
gpu::SparseContext<double> spmv_dev(ctx);   // share gpu::Context for same-stream
auto d_A = spmv_dev.deviceView(A);          // upload sparse matrix once
d_y = d_A * d_x;                            // operator syntax, stays on device
```

### Eigen algorithm interop (example: Conjugate gradient)

The BLAS-1 operators and `DeviceSparseView` make `DeviceMatrix` usable as a
vector type in GPU implementations of algorithms like conjugate gradient.
Conjugate gradient is the motivating example -- the GPU CG mirrors Eigen's
`conjugate_gradient()` line for line, with only one host sync per iteration
(the convergence check). All scalar intermediates (`alpha`, `beta`, `absNew`)
stay on device as `DeviceScalar` values:

```cpp
gpu::Context ctx;
gpu::Context::setThreadLocal(&ctx);
gpu::SparseContext<double> spmv(ctx);
auto mat = spmv.deviceView(A);              // upload sparse matrix once

auto rhs = gpu::DeviceMatrix<double>::fromHost(b, ctx.stream());
gpu::DeviceMatrix<double> x(n, 1);
x.setZero();
gpu::DeviceMatrix<double> residual(n, 1);
residual.copyFrom(ctx, rhs);                // r = b (x=0)
gpu::DeviceMatrix<double> p(n, 1);
p.copyFrom(ctx, residual);                  // p = r
gpu::DeviceMatrix<double> z(n, 1), tmp(n, 1);

auto absNew = residual.dot(p);              // DeviceScalar -- no sync

while (i < maxIters) {
  tmp.noalias() = mat * p;                   // SpMV, device-resident

  auto alpha = absNew / p.dot(tmp);          // DeviceScalar / DeviceScalar -- no sync

  x += alpha * p;                            // DeviceScalar * DeviceMatrix axpy -- no sync
  residual -= alpha * tmp;                   // DeviceScalar * DeviceMatrix axpy -- no sync

  residualNorm2 = residual.squaredNorm();    // THE one sync per iteration
  if (residualNorm2 < threshold) break;

  z.copyFrom(ctx, residual);                 // no preconditioner: z = r
  auto absOld = std::move(absNew);           // no sync, no alloc
  absNew = residual.dot(z);                  // DeviceScalar -- no sync
  auto beta = absNew / absOld;               // DeviceScalar / DeviceScalar -- no sync

  p *= beta;                                 // DeviceScalar scal -- no sync
  p += z;                                    // axpy -- no sync
}
MatrixXd result = x.toHost();
```

### Precision control

GEMM dispatch routes through `cublasLtMatmul`. The compute type is selected
per scalar via the `cuda_compute_type` trait in `CuBlasSupport.h`, gated by
two compile-time macros:

| Macro | Effect |
|---|---|
| (default) | `CUBLAS_COMPUTE_32F` / `CUBLAS_COMPUTE_64F`. cublasLt heuristics may pick tensor-core algorithms; on `sm_80+` doubles can land on Ozaki-emulated tensor cores. |
| `EIGEN_CUDA_TF32` | `CUBLAS_COMPUTE_32F_FAST_TF32` for `float` and `complex<float>` (~2x faster, 10-bit mantissa). No effect on `double` / `complex<double>`. |
| `EIGEN_NO_CUDA_TENSOR_OPS` | Pedantic compute types (`CUBLAS_COMPUTE_*_PEDANTIC`) for every scalar — disables tensor-core algorithms. Use for bit-exact reproducibility. Takes precedence over `EIGEN_CUDA_TF32`. |

These are independent of cuBLAS's runtime `cublasSetMathMode()` /
`CUBLAS_TF32_OVERRIDE` controls; the cublasLt path keys off the compile-time
compute type instead. The `cublasGemmEx` fallback (used when cublasLt's
heuristic returns no candidate) honors `EIGEN_NO_CUDA_TENSOR_OPS` via its
algorithm hint (`CUBLAS_GEMM_DEFAULT` vs `CUBLAS_GEMM_DEFAULT_TENSOR_OP`).

### Stream control and async execution

Operations are asynchronous by default. The compute-solve chain runs without
host synchronization until you need a result on the host:

```text
fromHost(A) --sync-->  compute() --async-->  solve() --async-->  toHost()
   H2D                  potrf                 potrs                D2H
                                                                   sync
```

Mandatory sync points:
- `fromHost()` -- Synchronizes to complete the upload before returning
- `toHost()` / `HostTransfer::get()` -- Must deliver data to host
- `info()` -- Must read the factorization status
- `DeviceScalar` implicit conversion -- Downloads scalar from device

**Cross-stream safety** is automatic. `DeviceMatrix` tracks write completion
via CUDA events. When a matrix written on stream A is read on stream B, the
module automatically inserts `cudaStreamWaitEvent`. Same-stream operations
skip the wait (CUDA guarantees in-order execution within a stream).

**Lifetime of cached factorizations.** A `gpu::LLT` / `gpu::LU` object owns
its CUDA stream, library handle, and the cached factor on device. Destroying
the factorization object while a `solve()` it issued is still in flight is
*correct* but not actually async: `cudaStreamDestroy` returns immediately,
but the destructor of the cached factor calls `cudaFree`, which is fully
device-synchronous and stalls until the in-flight `potrs`/`getrs` retires.
For genuine async pipelining keep the factorization object alive until you
have drained its results (e.g. via `toHost()` or by binding consumption to
an explicit `gpu::Context` that outlives both producer and consumer):

```cpp
gpu::LLT<double> llt(d_A);             // factor stays on device
auto h_x = d_X = llt.solve(d_B).toHostAsync(stream);
h_x.get();                             // sync: factor + result complete
// llt may now be destroyed without stalling the device
```

## Reference

### Supported scalar types

`float`, `double`, `std::complex<float>`, `std::complex<double>` (unless
noted otherwise).

### Expression -> library call mapping

| DeviceMatrix expression | Library call | Parameters |
|---|---|---|
| `C = A * B` | `cublasLtMatmul` (with `cublasGemmEx` fallback) | transA=N, transB=N, alpha=1, beta=0 |
| `C = A.adjoint() * B` | `cublasLtMatmul` | transA=C, transB=N |
| `C = A.transpose() * B` | `cublasLtMatmul` | transA=T, transB=N |
| `C = A * B.adjoint()` | `cublasLtMatmul` | transA=N, transB=C |
| `C = A * B.transpose()` | `cublasLtMatmul` | transA=N, transB=T |
| `C = alpha * A * B` | `cublasLtMatmul` | alpha from LHS |
| `C = A * (alpha * B)` | `cublasLtMatmul` | alpha from RHS |
| `C += A * B` | `cublasLtMatmul` | alpha=1, beta=1 |
| `C.device(ctx) -= A * B` | `cublasLtMatmul` | alpha=-1, beta=1 |
| `X = A.llt().solve(B)` | `cusolverDnXpotrf` + `Xpotrs` | uplo, n, nrhs |
| `X = A.llt<Upper>().solve(B)` | same | uplo=Upper |
| `X = A.lu().solve(B)` | `cusolverDnXgetrf` + `Xgetrs` | n, nrhs |
| `X = A.triangularView<L>().solve(B)` | `cublasXtrsm` | side=L, uplo, diag=NonUnit |
| `C = A.selfadjointView<L>() * B` | `cublasXsymm` / `cublasXhemm` | side=L, uplo |
| `C.selfadjointView<L>().rankUpdate(A)` | `cublasXsyrk` / `cublasXherk` | uplo, trans=N |
| `C = A + B` | `cublasXgeam` | alpha=1, beta=1 |
| `C = A + alpha * B` | `cublasXgeam` | alpha=1, beta from scaled |
| `C = A - B` | `cublasXgeam` | alpha=1, beta=-1 |
| `C = A - alpha * B` | `cublasXgeam` | alpha=1, beta=-scaled |
| `x += alpha * y` | `cublasXaxpy` | alpha (host scalar) |
| `x += dAlpha * y` | `cublasXaxpy` | alpha (DeviceScalar, device pointer mode) |
| `x -= alpha * y` | `cublasXaxpy` | alpha negated |
| `x *= alpha` | `cublasXscal` | alpha (host or DeviceScalar) |
| `x.dot(y)` | `cublasXdot` / `cublasXdotc` | returns `DeviceScalar` |
| `x.norm()` | `cublasXnrm2` | returns `DeviceScalar<RealScalar>` |
| `x.squaredNorm()` | `cublasXdot(x, x)` | returns `DeviceScalar<RealScalar>` |
| `d_y = view * d_x` | `cusparseSpMV` | device-resident SpMV |

### `DeviceMatrix<Scalar>`

Typed RAII wrapper for a dense column-major matrix in GPU device memory.
Always dense (leading dimension = rows). A vector is a `DeviceMatrix` with
one column.

```cpp
// Construction
DeviceMatrix<Scalar>()                                   // Empty (0x0)
DeviceMatrix<Scalar>(Index n)                            // Allocate column vector (n x 1)
DeviceMatrix<Scalar>(rows, cols)                         // Allocate uninitialized

// Upload / download / pointer adoption
static DeviceMatrix fromHost(matrix, stream=nullptr)           // -> DeviceMatrix (syncs)
static DeviceMatrix fromHostAsync(ptr, rows, cols, stream)         // -> DeviceMatrix (no sync, caller manages ptr lifetime)
static DeviceMatrix adopt(Scalar* device_ptr, rows, cols)          // Owning wrapper over a raw device pointer
static DeviceMatrix view(Scalar* device_ptr, rows, cols)           // Non-owning view (does not free on destruction)
PlainMatrix        toHost(stream=nullptr)                      // -> host Matrix (syncs)
HostTransfer       toHostAsync(stream=nullptr)                 // -> HostTransfer future (no sync)
DeviceMatrix       clone(stream=nullptr)                       // -> DeviceMatrix (D2D copy, async)

// Dimensions and access
Index   rows()
Index   cols()
size_t  sizeInBytes()
bool    empty()
Scalar* data()                                           // Raw device pointer
void    resize(Index rows, Index cols)                   // Discard contents, reallocate

// Expression builders (return lightweight views, evaluated on assignment)
AdjointView       adjoint()                              // GEMM with ConjTrans
TransposeView     transpose()                            // GEMM with Trans
LltExpr            llt() / llt<UpLo>()                   // -> .solve(d_B) -> DeviceMatrix
LuExpr             lu()                                  // -> .solve(d_B) -> DeviceMatrix
TriangularView     triangularView<UpLo>()                // -> .solve(d_B) -> DeviceMatrix (TRSM)
SelfAdjointView    selfadjointView<UpLo>()               // -> * d_B (SYMM), .rankUpdate(d_A) (SYRK)
Assignment   device(gpu::Context& ctx)                // Bind assignment to explicit stream
DeviceMatrix&      noalias()                             // No-op (all ops are implicitly noalias)

// BLAS Level-1 (all have overloads with explicit gpu::Context& parameter)
DeviceScalar<Scalar>     dot(const DeviceMatrix& other)  // cuBLAS dot/dotc -> DeviceScalar
DeviceScalar<RealScalar> norm()                          // cuBLAS nrm2 -> DeviceScalar
DeviceScalar<RealScalar>  squaredNorm()                    // dot(self, self) -> DeviceScalar (no sync)
void                     setZero()                       // cudaMemsetAsync
void                     addScaled(gpu::Context&, Scalar alpha, const DeviceMatrix& x)  // this += alpha * x (axpy)
void                     scale(gpu::Context&, Scalar alpha)                              // this *= alpha (scal)
void                     copyFrom(gpu::Context&, const DeviceMatrix& other)              // this = other (D2D copy)
DeviceMatrix& operator+=(const Scaled<DeviceMatrix>&)    // axpy; spelled `mat += alpha * other`
DeviceMatrix& operator-=(const Scaled<DeviceMatrix>&)    // axpy negated; spelled `mat -= alpha * other`
DeviceMatrix& operator+=(const DeviceMatrix&)            // cuBLAS axpy (alpha=1)
DeviceMatrix& operator-=(const DeviceMatrix&)            // cuBLAS axpy (alpha=-1)
DeviceMatrix& operator+=(const DeviceScaledDevice<Scalar>&)  // axpy with device scalar; spelled `mat += d_alpha * other`
DeviceMatrix& operator-=(const DeviceScaledDevice<Scalar>&)  // negated; spelled `mat -= d_alpha * other`
DeviceMatrix& operator*=(Scalar)                         // cuBLAS scal (host pointer mode)
DeviceMatrix& operator*=(const DeviceScalar<Scalar>&)    // cuBLAS scal (device pointer mode, no host sync)
DeviceMatrix  cwiseProduct(gpu::Context&, const DeviceMatrix&)            // NPP nppsMul (float/double only)
void          cwiseProduct(gpu::Context&, const DeviceMatrix&, const DeviceMatrix&)  // in-place: this = a .* b

// geam expressions (evaluated on assignment)
DeviceMatrix& operator=(const DeviceAddExpr&)            // C = A + B, C = A + alpha*B, C = A - B, etc.
```

### `DeviceScalar<Scalar>`

Device-resident scalar. Returned by `dot()`, `norm()`, and `squaredNorm()`.
Implicit conversion to `Scalar` triggers `cudaStreamSynchronize` + download.

```cpp
DeviceScalar(cudaStream_t stream = nullptr)              // Allocate uninitialized
DeviceScalar(Scalar host_val, cudaStream_t stream)       // Upload host value

Scalar         get()                                     // Download (syncs stream)
               operator Scalar()                         // Implicit conversion (syncs)
Scalar*        devicePtr()                               // Raw device pointer
cudaStream_t   stream()

// Device-side arithmetic (no host sync, real types only)
DeviceScalar   operator/(DeviceScalar, DeviceScalar)     // NPP nppsDiv
DeviceScalar   operator/(Scalar, DeviceScalar)           // upload + div
DeviceScalar   operator/(DeviceScalar, Scalar)           // upload + div
DeviceScalar   operator-()                               // NPP nppsMulC(-1)
```

### `gpu::Context`

Unified GPU execution context owning a CUDA stream and library handles. Not
thread-safe -- use one `Context` per thread, or external synchronization
across threads.

```cpp
gpu::Context()                                             // Creates dedicated stream + cuBLAS handle
                                                           // (cuSOLVER / cuBLASLt / cuSPARSE handles
                                                           // are created lazily on first use)
gpu::Context(cudaStream_t stream)                          // Borrow existing stream (not owned)
static gpu::Context& threadLocal()                         // Per-thread default (lazy-created)
static void        setThreadLocal(gpu::Context* ctx)       // Override thread-local default (nullptr restores)

cudaStream_t       stream()
cublasHandle_t     cublasHandle()
cusolverDnHandle_t cusolverHandle()                        // Lazy: creates the handle on first call
cublasLtHandle_t   cublasLtHandle()                        // Lazy-initialized
cusparseHandle_t   cusparseHandle()                        // Lazy-initialized

internal::DeviceBuffer*       gemmWorkspace()              // cublasLtMatmul scratch (lazy-grown per context)
internal::CublasLtPlanCache*  gemmPlanCache()              // shape-keyed plan cache (per context, ~8-entry LRU)
```

Non-copyable, non-movable (owns library handles). Translation units that
never call `cusolverHandle()` do not pull cuSOLVER symbols at link time --
see [Linking](#eigen_gpu_linking).

### `gpu::LLT<Scalar, UpLo>` -- Dense Cholesky (cuSOLVER)

Caches the Cholesky factor on device for repeated solves.

```cpp
gpu::LLT()                                                // Default construct, then call compute()
gpu::LLT(const EigenBase<D>& A)                           // Convenience: upload + factorize

gpu::LLT&            compute(const EigenBase<D>& A)       // Upload + factorize
gpu::LLT&            compute(const DeviceMatrix& d_A)     // D2D copy + factorize
gpu::LLT&            compute(DeviceMatrix&& d_A)          // Adopt + factorize (no copy)

PlainMatrix        solve(const MatrixBase<D>& B)         // -> host Matrix (syncs)
DeviceMatrix       solve(const DeviceMatrix& d_B)        // -> DeviceMatrix (async, stays on device)

ComputationInfo    info()                                // Lazy sync on first call: Success or NumericalIssue
Index              rows() / cols()
cudaStream_t       stream()
```

### `gpu::LU<Scalar>` -- Dense LU (cuSOLVER)

Same pattern as `gpu::LLT`. Adds a `gpu::GpuOp` parameter on `solve()`.

```cpp
PlainMatrix        solve(const MatrixBase<D>& B, GpuOp op = GpuOp::NoTrans)  // -> host Matrix
DeviceMatrix       solve(const DeviceMatrix& d_B, GpuOp op = GpuOp::NoTrans) // -> DeviceMatrix
```

`gpu::GpuOp`: `NoTrans`, `Trans`, `ConjTrans`.

### `gpu::QR<Scalar>` -- Dense QR (cuSOLVER)

QR factorization via `cusolverDnXgeqrf`. Solve uses ORMQR (apply Q^H) + TRSM
(back-substitute on R) -- Q is never formed explicitly.

```cpp
gpu::QR()                                                  // Default construct
gpu::QR(const EigenBase<D>& A)                             // Convenience: upload + factorize

gpu::QR&             compute(const EigenBase<D>& A)        // Upload + factorize
gpu::QR&             compute(const DeviceMatrix& d_A)      // D2D copy + factorize

PlainMatrix        solve(const MatrixBase<D>& B)         // -> host Matrix (syncs)
DeviceMatrix       solve(const DeviceMatrix& d_B)        // -> DeviceMatrix (async)
PlainMatrix        matrixR()                             // -> host Matrix (m >= n only)

ComputationInfo    info()                                // Lazy sync
Index              rows() / cols()
cudaStream_t       stream()
```

### `gpu::SVD<Scalar>` -- Dense SVD (cuSOLVER)

SVD via `cusolverDnXgesvd`. Supports `ComputeThinU | ComputeThinV`,
`ComputeFullU | ComputeFullV`, or `0` (values only). Wide matrices (m < n)
handled by internal transpose.

```cpp
gpu::SVD()                                                 // Default construct, then call compute()
gpu::SVD(const EigenBase<D>& A, unsigned options = ComputeThinU | ComputeThinV)  // Convenience

gpu::SVD&            compute(const EigenBase<D>& A, unsigned options = ComputeThinU | ComputeThinV)
gpu::SVD&            compute(const DeviceMatrix& d_A, unsigned options = ComputeThinU | ComputeThinV)

RealVector         singularValues()                      // -> host vector (syncs, downloads)
PlainMatrix        matrixU()                             // -> host Matrix (syncs, downloads)
PlainMatrix        matrixV()                             // -> host Matrix (V = VT^H, matches JacobiSVD)
PlainMatrix        matrixVT()                            // -> host Matrix (syncs, downloads V^T)

DeviceMatrix       d_singularValues()                    // -> DeviceMatrix view (zero-copy)
DeviceMatrix       d_matrixU()                           // -> DeviceMatrix view (zero-copy when m >= n)
DeviceMatrix       d_matrixVT()                          // -> DeviceMatrix view (zero-copy when m >= n)

PlainMatrix        solve(const MatrixBase<D>& B)         // -> host Matrix (pseudoinverse)
PlainMatrix        solve(const MatrixBase<D>& B, Index k)       // Truncated (top k triplets)
PlainMatrix        solve(const MatrixBase<D>& B, RealScalar l)  // Tikhonov regularized

Index              rank(RealScalar threshold = -1)
ComputationInfo    info()                                // Lazy sync
Index              rows() / cols()
cudaStream_t       stream()
```

**Note:** `singularValues()`, `matrixU()`, `matrixV()`, and `matrixVT()`
download to host on each call. The `d_*` accessors return non-owning
`DeviceMatrix` views into the solver's internal buffers; the `gpu::SVD` object
must outlive any view derived from it. For wide matrices (m < n) the U/V^T
views are owning (one `cublasXgeam` adjoint pass).

### `gpu::SelfAdjointEigenSolver<Scalar>` -- Eigendecomposition (cuSOLVER)

Symmetric/Hermitian eigenvalue decomposition via `cusolverDnXsyevd`.
`ComputeMode` enum: `EigenvaluesOnly`, `ComputeEigenvectors`.

```cpp
gpu::SelfAdjointEigenSolver()                              // Default construct, then call compute()
gpu::SelfAdjointEigenSolver(const EigenBase<D>& A, ComputeMode mode = ComputeEigenvectors)  // Convenience

gpu::SelfAdjointEigenSolver& compute(const EigenBase<D>& A, ComputeMode mode = ComputeEigenvectors)
gpu::SelfAdjointEigenSolver& compute(const DeviceMatrix& d_A, ComputeMode mode = ComputeEigenvectors)

RealVector         eigenvalues()                         // -> host vector (syncs, downloads, ascending order)
PlainMatrix        eigenvectors()                        // -> host Matrix (syncs, downloads, columns)

DeviceMatrix       d_eigenvalues()                       // -> DeviceMatrix view (zero-copy)
DeviceMatrix       d_eigenvectors()                      // -> DeviceMatrix view (zero-copy, requires ComputeEigenvectors)

ComputationInfo    info()                                // Lazy sync
Index              rows() / cols()
cudaStream_t       stream()
```

**Note:** `eigenvalues()` and `eigenvectors()` download to host on each call.
The `d_*` accessors return non-owning `DeviceMatrix` views into the solver's
internal buffers; the `gpu::SelfAdjointEigenSolver` object must outlive any
view derived from it.

### `HostTransfer<Scalar>`

Future for async device-to-host transfer. Returned by
`DeviceMatrix::toHostAsync()`.

```cpp
PlainMatrix&       get()                                 // Block until complete, return host Matrix ref. Idempotent.
bool               ready()                               // Non-blocking poll
```

### `gpu::SparseLLT<Scalar, UpLo>` -- Sparse Cholesky (cuDSS)

Requires cuDSS (CUDA 12.0+, `#define EIGEN_CUDSS`). Three-phase workflow
with symbolic reuse. Accepts `SparseMatrix<Scalar, ColMajor, int>` (CSC).
Matrix dimensions and nonzero count must fit in `int` (cuDSS limitation;
debug builds assert).

```cpp
gpu::SparseLLT()                                           // Default construct
gpu::SparseLLT(const SparseMatrixBase<D>& A)               // Analyze + factorize

gpu::SparseLLT&      analyzePattern(const SparseMatrixBase<D>& A)  // Symbolic analysis (reusable)
gpu::SparseLLT&      factorize(const SparseMatrixBase<D>& A)       // Numeric factorization
gpu::SparseLLT&      compute(const SparseMatrixBase<D>& A)         // analyzePattern + factorize

DenseMatrix        solve(const MatrixBase<D>& B)         // -> host Matrix (syncs)

ComputationInfo    info()                                // Lazy sync
Index              rows() / cols()
cudaStream_t       stream()
```

### `gpu::SparseLDLT<Scalar, UpLo>` -- Sparse LDL^T (cuDSS)

Symmetric indefinite. Same API as `gpu::SparseLLT`.

### `gpu::SparseLU<Scalar>` -- Sparse LU (cuDSS)

General non-symmetric. Same API as `gpu::SparseLLT` (without `UpLo`).

### `gpu::FFT<Scalar>` -- FFT (cuFFT)

Plans cached by (size, type) in a bounded LRU and reused; the
least-recently-used plan is destroyed via `cufftDestroy` on overflow. Cache
capacity is set at construction (default
`kDefaultCufftPlanCacheCapacity = 16`). Inverse transforms scaled so
`inv(fwd(x)) == x`. Supported scalars: `float`, `double`. Stream and cuBLAS
handle borrowed from a `gpu::Context` (default: `Context::threadLocal()`),
so by default the FFT shares a stream with other GPU operations on the same
thread.

```cpp
gpu::FFT(std::size_t plan_cache_capacity = kDefaultCufftPlanCacheCapacity)
                                        // bind to Context::threadLocal()
gpu::FFT(gpu::Context& ctx,
         std::size_t plan_cache_capacity = kDefaultCufftPlanCacheCapacity)
                                        // bind to an explicit Context

// 1D transforms (host vectors in and out)
ComplexVector      fwd(const MatrixBase<D>& x)           // C2C forward (complex input)
ComplexVector      fwd(const MatrixBase<D>& x)           // R2C forward (real input, returns n/2+1)
ComplexVector      inv(const MatrixBase<D>& X)           // C2C inverse, scaled by 1/n
RealVector         invReal(const MatrixBase<D>& X, Index n)  // C2R inverse, scaled by 1/n

// 2D transforms (host matrices in and out)
ComplexMatrix      fwd2(const MatrixBase<D>& A)         // 2D C2C forward
ComplexMatrix      inv2(const MatrixBase<D>& A)         // 2D C2C inverse, scaled by 1/(rows*cols)

cudaStream_t       stream()             // borrowed from the bound Context
gpu::Context&      context()            // the bound Context
std::size_t        plan_cache_capacity() // configured capacity
std::size_t        plan_cache_size()     // currently cached plan count
```

All FFT methods accept host data and return host data. Upload/download is
handled internally. The C2C and R2C overloads of `fwd()` are distinguished by
the input scalar type (complex vs real).

### `gpu::SparseContext<Scalar>` -- SpMV/SpMM (cuSPARSE)

Accepts `SparseMatrix<Scalar, ColMajor>`. Host-input methods accept host data
and return host data; device-input methods (`deviceView()`, `multiply(A, d_x,
d_y)`) operate on `DeviceMatrix`. Matrix dimensions and nonzero count must fit
in `int` (cuSPARSE limitation; debug builds assert).

```cpp
gpu::SparseContext()                                       // Creates own stream + cuSPARSE handle
gpu::SparseContext(gpu::Context& ctx)                        // Borrow gpu::Context for same-stream execution

// Host data in/out
DenseVector        multiply(A, x)                                       // y = A * x
void               multiply(A, x, y, alpha=1, beta=0,                   // y = alpha*op(A)*x + beta*y
                     op=GpuOp::NoTrans)
DenseVector        multiplyT(A, x)                                      // y = A^T * x
DenseVector        multiplyAdjoint(A, x)                                // y = A^H * x
DenseMatrix        multiplyMat(A, X, op=GpuOp::NoTrans)                 // Y = op(A) * X (SpMM)

// DeviceMatrix in/out (sparse matrix re-uploaded each call)
void               multiply(A, d_x, d_y)                                // SpMV with device vectors
void               multiply(A, d_x, d_y, alpha, beta, op)

// Device-resident sparse matrix (upload once, reuse)
DeviceSparseView   deviceView(A)                                        // Upload sparse matrix, return view

cudaStream_t       stream()
```

### `DeviceSparseView<Scalar>` -- Device-resident sparse matrix

Returned by `gpu::SparseContext::deviceView()`. Holds a sparse matrix on device
for repeated SpMV without re-uploading.

```cpp
SpMVExpr           operator*(const DeviceMatrix& d_x)    // d_y = view * d_x (evaluated on assignment)
```

### Aliasing

Unlike Eigen's `Matrix`, where omitting `.noalias()` triggers a copy to a
temporary, DeviceMatrix dispatches directly to NVIDIA library calls which have
no built-in aliasing protection. All operations are implicitly noalias.
The caller must ensure operands don't alias the destination for GEMM, TRSM,
SYMM/HEMM, and SYRK/HERK. Debug builds assert on these violations before
dispatching to cuBLAS. `geam` expressions (`d_C = d_A + alpha * d_B`) are
safe with aliasing. The `.noalias()` method exists as a no-op for Eigen
template compatibility.

## Future work

- **Reassess host-input vs. device-input API surface.** Each solver currently
  exposes both host-input (`compute(MatrixXd)`, `solve(MatrixXd)`) and
  device-input (`compute(DeviceMatrix)`, `solve(DeviceMatrix)`) overloads, plus
  host- and device-side accessors (`matrixU()` vs `d_matrixU()`). This eases
  migration from CPU Eigen but may invite accidental host ↔ device round-trips
  when users mix the two without realising the cost. Revisit once the module
  is in users' hands; if the convenience overloads cause more confusion than
  they save, narrow toward a single explicit `fromHost` / `toHost` boundary.

- **cuDSS configuration knobs.** cuDSS exposes settings for accuracy /
  robustness (e.g. matching, pivoting) and execution mode (e.g. hybrid
  memory, hybrid execute). The current bindings use cuDSS defaults, which
  are tuned for performance rather than maximum robustness — for example,
  matching is off by default. We don't expose configuration controls yet;
  a follow-up should add a `gpu::SparseSolverConfig` (or per-solver
  setters) covering at least matching, pivot threshold, and reordering
  algorithm pass-through, and consider switching the defaults toward
  robustness once exposed.

- **cuDSS threading layer for host-side reordering.** As of cuDSS 0.7.1
  fill-reducing reordering runs on the CPU. cuDSS supports a "threading
  layer" plugin that parallelises this stage; for reordering-dominated
  problems it can materially close the gap with multithreaded CPU sparse
  direct solvers. We don't currently configure a threading layer.

- **Complex symmetric (non-Hermitian) sparse LDL^T.** `gpu::SparseLDLT`
  treats complex inputs as Hermitian (matching `Eigen::SimplicialLDLT`).
  cuDSS also supports `CUDSS_MTYPE_SYMMETRIC` for complex matrices
  (A = A^T, no conjugation); exposing this would need a separate solver
  mode.

## File layout

| File | Depends on | Contents |
|------|-----------|----------|
| `GpuSupport.h` | `<cuda_runtime.h>` | Error macro, `DeviceBuffer`, `DeviceBufferPool`, `cuda_data_type<>` |
| `DeviceMatrix.h` | `GpuSupport.h` | `gpu::DeviceMatrix<>`, `gpu::HostTransfer<>` |
| `DeviceExpr.h` | `DeviceMatrix.h` | GEMM, geam, and device-scalar expression wrappers |
| `DeviceBlasExpr.h` | `DeviceMatrix.h` | TRSM, SYMM, SYRK expression wrappers |
| `DeviceSolverExpr.h` | `DeviceMatrix.h` | Solver expression wrappers (LLT, LU) |
| `DeviceScalar.h` | `GpuSupport.h`, `DeviceScalarOps.h` | `gpu::DeviceScalar<>` (device-resident scalar) |
| `DeviceScalarOps.h` | `<npps_*.h>` | Scalar div/neg/cwiseProduct via NPP |
| `DeviceDispatch.h` | all above | All dispatch functions, BLAS-1 out-of-line defs, `gpu::Assignment` |
| `GpuContext.h` | `CuBlasSupport.h`, `CuSolverSupport.h` | `gpu::Context` |
| `CuBlasSupport.h` | `GpuSupport.h`, `<cublas_v2.h>`, `<cublasLt.h>` | cuBLAS error macro, type-specific wrappers |
| `CuSolverSupport.h` | `GpuSupport.h`, `<cusolverDn.h>` | cuSOLVER params, fill-mode mapping |
| `GpuSolverContext.h` | `CuSolverSupport.h`, `CuBlasSupport.h` | Shared solver context (stream, handles, scratch) |
| `GpuLLT.h` | `GpuSolverContext.h` | `gpu::LLT<>` -- Cached dense Cholesky factorization |
| `GpuLU.h` | `GpuSolverContext.h` | `gpu::LU<>` -- Cached dense LU factorization |
| `GpuQR.h` | `GpuSolverContext.h` | `gpu::QR<>` -- Dense QR decomposition |
| `GpuSVD.h` | `GpuSolverContext.h` | `gpu::SVD<>` -- Dense SVD decomposition |
| `GpuEigenSolver.h` | `GpuSolverContext.h` | `gpu::SelfAdjointEigenSolver<>` |
| `CuFftSupport.h` | `GpuSupport.h`, `<cufft.h>` | cuFFT error macro, type-dispatch wrappers |
| `GpuFFT.h` | `CuFftSupport.h`, `CuBlasSupport.h`, `GpuContext.h` | `gpu::FFT<>` -- 1D/2D FFT with plan caching |
| `CuSparseSupport.h` | `GpuSupport.h`, `<cusparse.h>` | cuSPARSE error macro |
| `GpuSparseContext.h` | `CuSparseSupport.h` | `gpu::SparseContext<>`, `gpu::DeviceSparseView<>` |
| `CuDssSupport.h` | `GpuSupport.h`, `<cudss.h>` | cuDSS error macro, type traits (optional) |
| `GpuSparseSolverBase.h` | `CuDssSupport.h` | CRTP base for sparse solvers (optional) |
| `GpuSparseLLT.h` | `GpuSparseSolverBase.h` | `gpu::SparseLLT<>` -- Sparse Cholesky via cuDSS (optional) |
| `GpuSparseLDLT.h` | `GpuSparseSolverBase.h` | `gpu::SparseLDLT<>` -- Sparse LDL^T via cuDSS (optional) |
| `GpuSparseLU.h` | `GpuSparseSolverBase.h` | `gpu::SparseLU<>` -- Sparse LU via cuDSS (optional) |

## Building and testing

```bash
cmake -G Ninja -B build -S . \
  -DEIGEN_TEST_CUDA=ON \
  -DEIGEN_CUDA_COMPUTE_ARCH="70" \
  -DEIGEN_TEST_CUBLAS=ON \
  -DEIGEN_TEST_CUSOLVER=ON

cmake --build build --target cublas cusolver_llt cusolver_lu \
  cusolver_qr cusolver_svd cusolver_eigen \
  device_matrix cufft cusparse_spmv cg
ctest --test-dir build -L gpu --output-on-failure

# Sparse solvers (cuDSS -- separate install required)
cmake -G Ninja -B build -S . \
  -DEIGEN_TEST_CUDA=ON \
  -DEIGEN_CUDA_COMPUTE_ARCH="70" \
  -DEIGEN_TEST_CUDSS=ON

cmake --build build --target cudss_llt cudss_ldlt cudss_lu
ctest --test-dir build -R '^cudss_' --output-on-failure
```

## Future enhancements

- **Batched API (`DeviceBatchMatrix`).** A strided batch of N identical-size
  matrices dispatching to cuBLAS/cuSOLVER batched APIs (`cublasDgemmBatched`,
  `cusolverDnXpotrfBatched`, etc.). This enables robotics and model-predictive
  control workloads where many small independent systems are solved in
  parallel.
- **cuTENSOR for Tensor module.** Replace the hand-written GPU tensor
  contraction and reduction kernels (~2300 lines in
  `TensorContractionGpu.h` / `TensorReductionGpu.h`) with cuTENSOR dispatch,
  following the same library-dispatch pattern used by `unsupported/Eigen/GPU`.
- **Unified/zero-copy memory for Jetson.** Use `cudaMallocManaged` or
  `cudaHostAllocMapped` to eliminate `fromHost()` / `toHost()` copies on
  integrated GPUs (Jetson) where CPU and GPU share DRAM.
- **Device-side Eigen interop.** Bridge between host-side `DeviceMatrix`
  dispatch and device-side Eigen expression templates (Core + Tensor) running
  inside CUDA kernels. Raw-pointer + `Map` / `TensorMap` as the zero-copy
  interop surface.
- **Per-stream CUDA memory pools.** Currently all streams in a `GpuContext`
  share the default device memory pool. Attaching a dedicated
  `cudaMemPool_t` per stream (`cudaDeviceSetMempool` /
  `cudaMallocFromPoolAsync`) can reduce cross-stream allocator contention for
  workloads that fan out many concurrent solves.
