<!--
SPDX-FileCopyrightText: The Eigen Authors
SPDX-License-Identifier: MPL-2.0
-->

# Bundle Adjustment: GPU CG vs CPU CG Results

Benchmark of Eigen's GPU CG pipeline on normal equations arising from bundle
adjustment (BAL datasets). Compares CPU `ConjugateGradient` (Jacobi preconditioner)
against GPU CG using `DeviceMatrix` + `GpuSparseContext` + `DeviceScalar`.

## Hardware

- **CPU**: Intel Core i7-13700HX (Raptor Lake, 12 cores / 24 threads, single thread for Eigen CG)
- **GPU**: NVIDIA GeForce RTX 4070 Laptop GPU (Ada Lovelace, 4608 CUDA cores, 8 GB GDDR6)
- **CUDA**: 13.2 / Driver 595.79
- **OS**: Ubuntu 24.04 (WSL2, kernel 6.6.87)

## Software

- Eigen: `eigen-gpu-cg` branch
- Google Benchmark 1.9.1
- Compiler: nvcc 13.2 + g++ 13.3
- Normal equations: H = J^T*J + I (Levenberg-Marquardt damping lambda=1.0)
- CG tolerance: 1e-8, max iterations: 10000

## Method

For each BAL problem file:
1. Parse the BAL file (cameras, 3D points, 2D observations)
2. Compute the full Jacobian J using the BAL camera model (Rodrigues rotation +
   perspective projection + radial distortion) with central finite differences
3. Form the normal equations H = J^T*J + lambda*I (sparse, symmetric positive definite)
4. Solve H*dx = -J^T*r using CG with Jacobi preconditioner on CPU and GPU
5. Report wall-clock time (mean of 3 repetitions)

GPU CG uses: `GpuSparseContext` for SpMV, `DeviceMatrix` for vectors,
`DeviceScalar` with `CUBLAS_POINTER_MODE_DEVICE` for dot/norm reductions,
in-place `cwiseProduct` via NPP for Jacobi preconditioner application,
device-pointer-mode `scal` to avoid host sync on the beta update.

## Results

### Summary table

| Dataset | Cameras | Points | Obs | H size | H nnz | CG iters | CPU CG (ms) | GPU CG (ms) | Speedup |
|---------|---------|--------|-----|--------|-------|----------|-------------|-------------|---------|
| Ladybug-49 | 49 | 7,776 | 31,843 | 23,769 | 1.8M | 4,421 | 4,006 | 1,152 | **3.5x** |
| Ladybug-138 | 138 | 19,878 | 85,217 | 60,876 | 4.8M | 7,008 | 21,498 | 3,553 | **6.1x** |
| Ladybug-646 | 646 | 73,584 | 327,297 | 226,566 | 18.4M | 10,000* | 123,727 | 14,268 | **8.7x** |
| Dubrovnik-356 | 356 | 226,730 | 1,255,268 | 683,394 | 69.8M | 4,308 | 216,149 | 24,493 | **8.8x** |

\* Hit 10,000 iteration cap (poorly conditioned problem). Both CPU and GPU
hit the same cap, so timing comparison remains valid.

### Profile breakdown (Ladybug-138, nsys)

GPU kernel time is dominated by SpMV (91%). The remaining 9% is BLAS-1
operations (dot, axpy, scal) and NPP element-wise ops (cwiseProduct).

| Kernel | Time (ms) | % | Calls |
|--------|-----------|---|-------|
| cuSPARSE csrmv (SpMV) | 2507 | 91.3% | 7,006 |
| cuBLAS dot | 92 | 3.4% | 21,020 |
| cuBLAS axpy (device ptr) | 27 | 1.0% | 14,012 |
| cuSPARSE partition | 19 | 0.7% | 7,006 |
| NPP cwiseProduct | 16 + 13 | 1.1% | 14,011 + 7,006 |
| cuBLAS axpy (host ptr) | 12 | 0.5% | 7,005 |
| cuBLAS scal (device ptr) | 11 | 0.4% | 7,005 |
| NPP scalar ops | 7 | 0.2% | 7,006 |

### Optimizations applied

Three profiling-driven optimizations reduced GPU CG time by **1.8x**
(6.5s → 3.6s on Ladybug-138):

1. **In-place `cwiseProduct`**: The Jacobi preconditioner apply
   (`z = invdiag .* residual`) was allocating a new DeviceMatrix every
   iteration. Added `z.cwiseProduct(ctx, a, b)` that reuses `z`'s buffer.
   Reduced `cudaMalloc` calls from 7,053 to 23 (saving 2.3s).

2. **`squaredNorm` via `dot(x,x)`**: cuBLAS `nrm2` uses a numerically
   careful scaled-sum-of-squares algorithm (29µs/call). Replaced with
   `dot(x,x)` (6.4µs/call) — 4.5x faster per call, saving ~320ms.

3. **Device-pointer `scal`**: `p *= beta` was converting `DeviceScalar`
   beta to host (triggering a stream sync), then calling host-pointer-mode
   scal. Added `operator*=(DeviceScalar)` that uses device-pointer-mode
   scal, eliminating one sync per iteration. Halved `cudaStreamSynchronize`
   calls from 14K to 7K.

### Observations

1. **GPU speedup scales with problem size**: from 3.5x on small problems
   (24K variables) to 8.8x on large problems (683K variables). This is
   expected — larger problems have more parallelism for the GPU to exploit.

2. **Iteration counts match**: CPU and GPU CG converge in the same number
   of iterations (within 1%), confirming numerical equivalence.

3. **Bottleneck is SpMV**: CG iteration time is dominated (91%) by the
   sparse matrix-vector product on H. Further speedup requires either
   faster SpMV (e.g., block-sparse formats) or algorithmic improvements
   (Schur complement, better preconditioners).

4. **Remaining overhead**: CUDA API calls (cudaMemcpyAsync for 8-byte
   DeviceScalar transfers) account for ~50% of non-kernel time. Batching
   multiple scalar reductions into a single transfer would help.

5. **Jacobi preconditioner is weak for BA**: The Ladybug-646 problem does
   not converge in 10K iterations. Ceres uses block Jacobi or Schur
   complement preconditioners that would also benefit from GPU acceleration.

### Scaling plot data

```text
# n        nnz_H       cpu_ms    gpu_ms    speedup
23769      1793475     4006      1152      3.48
60876      4791762     21498     3553      6.05
226566     18387948    123727    14268     8.67
683394     69827066    216149    24493     8.82
```

## BAL datasets

Downloaded from http://grail.cs.washington.edu/projects/bal/

| File | Source |
|------|--------|
| problem-49-7776-pre.txt | Ladybug sequence |
| problem-138-19878-pre.txt | Ladybug sequence |
| problem-646-73584-pre.txt | Ladybug sequence |
| problem-356-226730-pre.txt | Dubrovnik reconstruction |

## Reproducing

```bash
# Build
cmake -G Ninja -B build-bench-gpu -S unsupported/benchmarks/GPU -DCMAKE_CUDA_ARCHITECTURES=89
cmake --build build-bench-gpu --target bench_ba

# Download BAL datasets
wget http://grail.cs.washington.edu/projects/bal/data/ladybug/problem-49-7776-pre.txt.bz2
wget http://grail.cs.washington.edu/projects/bal/data/ladybug/problem-138-19878-pre.txt.bz2
wget http://grail.cs.washington.edu/projects/bal/data/ladybug/problem-646-73584-pre.txt.bz2
wget http://grail.cs.washington.edu/projects/bal/data/dubrovnik/problem-356-226730-pre.txt.bz2
bunzip2 *.bz2

# Run (one at a time)
BAL_FILE=problem-49-7776-pre.txt ./build-bench-gpu/bench_ba --benchmark_repetitions=3
BAL_FILE=problem-138-19878-pre.txt ./build-bench-gpu/bench_ba --benchmark_repetitions=3
BAL_FILE=problem-646-73584-pre.txt ./build-bench-gpu/bench_ba --benchmark_repetitions=3
BAL_FILE=problem-356-226730-pre.txt ./build-bench-gpu/bench_ba --benchmark_repetitions=3
```
