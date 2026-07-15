# GPU benchmarks require CUDA runtime + cuSOLVER.
# Build separately from the main benchmark tree since they need CUDA toolchain.
#
# Usage:
#   cmake -G Ninja -B build-bench-gpu -S unsupported/benchmarks/GPU \
#         -DCMAKE_CUDA_ARCHITECTURES=89
#   cmake --build build-bench-gpu
#
# Profiling:
#   nsys profile --trace=cuda ./build-bench-gpu/bench_solvers
#   ncu --set full -o profile ./build-bench-gpu/bench_solvers --benchmark_filter=BM_GpuLLT_Compute/4096
# SPDX-FileCopyrightText: The Eigen Authors
# SPDX-License-Identifier: MPL-2.0

cmake_minimum_required(VERSION 3.18)
project(EigenGpuBenchmarks CXX CUDA)

find_package(benchmark REQUIRED)
find_package(CUDAToolkit REQUIRED)

set(EIGEN_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../..")

function(eigen_add_gpu_benchmark name source)
  cmake_parse_arguments(BENCH "" "" "LIBRARIES;DEFINITIONS" ${ARGN})
  if(NOT IS_ABSOLUTE "${source}")
    set(source "${CMAKE_CURRENT_SOURCE_DIR}/${source}")
  endif()
  add_executable(${name} ${source})
  target_compile_features(${name} PRIVATE cxx_std_14)
  target_include_directories(${name} PRIVATE
    ${EIGEN_SOURCE_DIR}
    ${CUDAToolkit_INCLUDE_DIRS})
  target_link_libraries(${name} PRIVATE
    benchmark::benchmark benchmark::benchmark_main
    CUDA::cudart CUDA::cusolver CUDA::cublas)
  if(BENCH_LIBRARIES)
    target_link_libraries(${name} PRIVATE ${BENCH_LIBRARIES})
  endif()
  # Note: do NOT set -DNDEBUG here. EIGEN_CUDA_RUNTIME_CHECK uses eigen_assert,
  # which is a no-op under NDEBUG, silently skipping cudaStreamSynchronize calls.
  target_compile_options(${name} PRIVATE -O3)
  target_compile_definitions(${name} PRIVATE EIGEN_USE_GPU)
  if(BENCH_DEFINITIONS)
    target_compile_definitions(${name} PRIVATE ${BENCH_DEFINITIONS})
  endif()
endfunction()

# Solver benchmarks: LLT/LU compute + solve, host vs device paths, CPU baselines.
eigen_add_gpu_benchmark(bench_solvers bench_solvers.cpp)
eigen_add_gpu_benchmark(bench_solvers_float bench_solvers.cpp DEFINITIONS SCALAR=float)

# GEMM benchmarks: cublasLtMatmul (with plan cache) vs raw cublasGemmEx across sizes.
eigen_add_gpu_benchmark(bench_gemm bench_gemm.cpp LIBRARIES CUDA::cublasLt)
eigen_add_gpu_benchmark(bench_gemm_float bench_gemm.cpp LIBRARIES CUDA::cublasLt DEFINITIONS SCALAR=float)

# Chaining benchmarks: async pipeline efficiency, host-roundtrip vs device chain.
eigen_add_gpu_benchmark(bench_chaining bench_chaining.cpp)
eigen_add_gpu_benchmark(bench_chaining_float bench_chaining.cpp DEFINITIONS SCALAR=float)

# Batching benchmarks: multi-stream concurrency for many small systems.
eigen_add_gpu_benchmark(bench_batching bench_batching.cpp)
eigen_add_gpu_benchmark(bench_batching_float bench_batching.cpp DEFINITIONS SCALAR=float)

# FFT benchmarks: 1D/2D C2C, R2C, C2R throughput and plan reuse.
eigen_add_gpu_benchmark(bench_fft bench_fft.cpp LIBRARIES CUDA::cufft)
eigen_add_gpu_benchmark(bench_fft_double bench_fft.cpp LIBRARIES CUDA::cufft DEFINITIONS SCALAR=double)

# CG sync overhead benchmark: host vs device pointer mode for reductions.
# Uses CUDA kernels for device scalar arithmetic.
add_executable(bench_cg_sync bench_cg_sync.cu)
target_include_directories(bench_cg_sync PRIVATE
  ${EIGEN_SOURCE_DIR}
  ${CUDAToolkit_INCLUDE_DIRS})
target_link_libraries(bench_cg_sync PRIVATE
  benchmark::benchmark benchmark::benchmark_main
  CUDA::cudart CUDA::cusolver CUDA::cublas CUDA::cusparse CUDA::npps CUDA::nppc)
target_compile_options(bench_cg_sync PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-O3 --expt-relaxed-constexpr>)
target_compile_definitions(bench_cg_sync PRIVATE EIGEN_USE_GPU)

# GPU CG vs CPU CG comparison benchmark.
add_executable(bench_cg_vs_cpu bench_cg_vs_cpu.cu)
target_include_directories(bench_cg_vs_cpu PRIVATE
  ${EIGEN_SOURCE_DIR}
  ${CUDAToolkit_INCLUDE_DIRS})
target_link_libraries(bench_cg_vs_cpu PRIVATE
  benchmark::benchmark benchmark::benchmark_main
  CUDA::cudart CUDA::cusolver CUDA::cublas CUDA::cusparse CUDA::npps CUDA::nppc)
target_compile_options(bench_cg_vs_cpu PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-O3 --expt-relaxed-constexpr>)
target_compile_definitions(bench_cg_vs_cpu PRIVATE EIGEN_USE_GPU)

# Bundle Adjustment benchmark: GPU CG vs CPU CG on real BAL datasets.
add_executable(bench_ba bench_ba.cu)
target_include_directories(bench_ba PRIVATE
  ${EIGEN_SOURCE_DIR}
  ${CUDAToolkit_INCLUDE_DIRS})
target_link_libraries(bench_ba PRIVATE
  benchmark::benchmark
  CUDA::cudart CUDA::cusolver CUDA::cublas CUDA::cusparse CUDA::npps CUDA::nppc)
target_compile_options(bench_ba PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-O3 --expt-relaxed-constexpr>)
target_compile_definitions(bench_ba PRIVATE EIGEN_USE_GPU)
