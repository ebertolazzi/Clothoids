/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2013-2025                                                 |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#pragma once

#ifndef SYSTEM_UTILS_HXX
#define SYSTEM_UTILS_HXX

#include "Utils.hh"

// Include for macOS specific APIs
#if defined( __APPLE__ )
#include <mach/mach.h>
#include <mach/mach_host.h>
#include <sys/sysctl.h>
#include <sys/mount.h>
#include <sys/statvfs.h>
#endif

namespace Utils
{

  /*==========================================================================*\
  |                            CLASS Architecture                              |
  \*==========================================================================*/

  /**
   * @class Architecture
   * @brief Provides architecture detection and CPU feature querying utilities.
   *
   * This class encapsulates functionality to detect the system architecture
   * (x86, ARM, PowerPC, etc.) and query CPU-specific features (MMX, SSE, AVX,
   * NEON, etc.) for performance optimization and compatibility checking.
   *
   * @note The class uses compile-time macros and runtime detection where needed.
   */
  class Architecture
  {
  public:
    /**
     * @enum ArchType
     * @brief Enumeration of supported CPU architectures.
     */
    enum class ArchType
    {
      UNKNOWN,    ///< Unknown architecture
      X86,        ///< 32-bit x86 architecture
      X64,        ///< 64-bit x86-64 architecture
      ARM32,      ///< 32-bit ARM architecture
      ARM64,      ///< 64-bit ARM architecture (AArch64)
      POWERPC,    ///< 32-bit PowerPC architecture
      POWERPC64,  ///< 64-bit PowerPC architecture
      MIPS,       ///< 32-bit MIPS architecture
      MIPS64,     ///< 64-bit MIPS architecture
      RISCV32,    ///< 32-bit RISC-V architecture
      RISCV64     ///< 64-bit RISC-V architecture
    };

    // =========================================================================
    // BASIC ARCHITECTURE DETECTION
    // =========================================================================

    /**
     * @brief Detects the current system architecture.
     * @return ArchType enumeration value representing the detected architecture.
     *
     * This function uses compiler-defined macros to determine the architecture
     * at compile time. It's useful for conditional compilation and runtime
     * feature selection.
     */
    static ArchType get_architecture()
    {
#if defined( __x86_64__ ) || defined( _M_X64 )
      return ArchType::X64;
#elif defined( __i386__ ) || defined( _M_IX86 )
      return ArchType::X86;
#elif defined( __aarch64__ ) || defined( _M_ARM64 )
      return ArchType::ARM64;
#elif defined( __arm__ ) || defined( _M_ARM )
      return ArchType::ARM32;
#elif defined( __powerpc64__ ) || defined( __ppc64__ )
      return ArchType::POWERPC64;
#elif defined( __powerpc__ ) || defined( __ppc__ )
      return ArchType::POWERPC;
#elif defined( __mips64__ )
      return ArchType::MIPS64;
#elif defined( __mips__ )
      return ArchType::MIPS;
#elif defined( __riscv ) && ( __riscv_xlen == 64 )
      return ArchType::RISCV64;
#elif defined( __riscv ) && ( __riscv_xlen == 32 )
      return ArchType::RISCV32;
#else
      return ArchType::UNKNOWN;
#endif
    }

    /**
     * @brief Returns a human-readable string for the current architecture.
     * @return String representation of the architecture.
     */
    static std::string get_architecture_string()
    {
      switch ( get_architecture() )
      {
        case ArchType::X64: return "x86_64";
        case ArchType::X86: return "x86";
        case ArchType::ARM64: return "ARM64";
        case ArchType::ARM32: return "ARM32";
        case ArchType::POWERPC64: return "PowerPC64";
        case ArchType::POWERPC: return "PowerPC";
        case ArchType::MIPS64: return "MIPS64";
        case ArchType::MIPS: return "MIPS";
        case ArchType::RISCV64: return "RISC-V64";
        case ArchType::RISCV32: return "RISC-V32";
        default: return "Unknown";
      }
    }

    // =========================================================================
    // X86/X64 CPU DETECTION
    // =========================================================================

#ifdef UTILS_CPU_X86

    // Bit definitions for CPUID feature flags in EDX register
    static constexpr unsigned long EDX_MMX_bit      = 0x800000UL;    ///< MMX technology (bit 23)
    static constexpr unsigned long EDX_SSE_bit      = 0x2000000UL;   ///< SSE instructions (bit 25)
    static constexpr unsigned long EDX_SSE2_bit     = 0x4000000UL;   ///< SSE2 instructions (bit 26)
    static constexpr unsigned long EDX_3DnowExt_bit = 0x40000000UL;  ///< Extended 3DNow! (bit 30)
    static constexpr unsigned long EDX_3Dnow_bit    = 0x80000000UL;  ///< 3DNow! instructions (bit 31)
    static constexpr unsigned long EDX_MMXplus_bit  = 0x400000UL;    ///< MMX+ extensions (bit 22)

    // Bit definitions for CPUID feature flags in ECX register
    static constexpr unsigned long ECX_SSE3_bit    = 0x1UL;         ///< SSE3 instructions (bit 0)
    static constexpr unsigned long ECX_SSSE3_bit   = 0x200UL;       ///< SSSE3 instructions (bit 9)
    static constexpr unsigned long ECX_SSE41_bit   = 0x80000UL;     ///< SSE4.1 instructions (bit 19)
    static constexpr unsigned long ECX_SSE42_bit   = 0x100000UL;    ///< SSE4.2 instructions (bit 20)
    static constexpr unsigned long ECX_SSE4A_bit   = 0x40UL;        ///< SSE4a instructions (AMD, bit 6)
    static constexpr unsigned long ECX_SSE5_bit    = 0x800UL;       ///< SSE5 instructions (AMD, bit 11)
    static constexpr unsigned long ECX_AVX_bit     = 0x10000000UL;  ///< AVX instructions (bit 28)
    static constexpr unsigned long ECX_AVX2_bit    = 0x20UL;        ///< AVX2 instructions (bit 5, CPUID level 7)
    static constexpr unsigned long ECX_AVX512F_bit = 0x10000UL;     ///< AVX-512 Foundation (bit 16, CPUID level 7)

    /**
     * @struct CPUIDResult
     * @brief Container for CPUID instruction results.
     */
    struct CPUIDResult
    {
      unsigned long eax;  ///< EAX register value
      unsigned long ebx;  ///< EBX register value
      unsigned long ecx;  ///< ECX register value
      unsigned long edx;  ///< EDX register value
    };

    /**
     * @brief Executes the CPUID instruction with the specified level.
     * @param[out] result Structure to store the CPUID results.
     * @param[in] level CPUID function/level to execute.
     */
    static void cpuId( CPUIDResult & result, unsigned long level )
    {
#if defined( _WIN32 ) || defined( _WIN64 )
      int cpuInfo[4];
      __cpuidex( cpuInfo, (int) level, 0 );
      result.eax = cpuInfo[0];
      result.ebx = cpuInfo[1];
      result.ecx = cpuInfo[2];
      result.edx = cpuInfo[3];
#else
#if defined( __i386__ ) && defined( __PIC__ )
      // Handle PIC register in ebx for 32-bit PIC code
      __asm__ volatile(
        "xchgl %%ebx, %1\n\t"
        "cpuid\n\t"
        "xchgl %%ebx, %1\n\t"
        : "=a"( result.eax ), "=r"( result.ebx ), "=c"( result.ecx ), "=d"( result.edx )
        : "0"( level ) );
#else
      __asm__ volatile( "cpuid\n\t"
                        : "=a"( result.eax ), "=b"( result.ebx ), "=c"( result.ecx ), "=d"( result.edx )
                        : "0"( level ) );
#endif
#endif
    }

    /**
     * @brief Executes the CPUID instruction with level and sublevel.
     * @param[out] result Structure to store the CPUID results.
     * @param[in] level CPUID function/level to execute.
     * @param[in] sublevel CPUID subfunction (in ECX register).
     */
    static void cpuId( CPUIDResult & result, unsigned long level, unsigned long sublevel )
    {
#if defined( _WIN32 ) || defined( _WIN64 )
      int cpuInfo[4];
      __cpuidex( cpuInfo, (int) level, (int) sublevel );
      result.eax = cpuInfo[0];
      result.ebx = cpuInfo[1];
      result.ecx = cpuInfo[2];
      result.edx = cpuInfo[3];
#else
      __asm__ volatile( "cpuid\n\t"
                        : "=a"( result.eax ), "=b"( result.ebx ), "=c"( result.ecx ), "=d"( result.edx )
                        : "0"( level ), "2"( sublevel ) );
#endif
    }

    /**
     * @struct CPUFeatures
     * @brief Holds x86/x64 CPU feature flags.
     */
    struct CPUFeatures
    {
      bool mmx           = false;  ///< MMX technology support
      bool mmxplus       = false;  ///< MMX+ extensions (AMD)
      bool sse           = false;  ///< SSE instructions
      bool sse2          = false;  ///< SSE2 instructions
      bool sse3          = false;  ///< SSE3 instructions
      bool ssse3         = false;  ///< SSSE3 instructions
      bool sse41         = false;  ///< SSE4.1 instructions
      bool sse42         = false;  ///< SSE4.2 instructions
      bool sse4a         = false;  ///< SSE4a instructions (AMD)
      bool sse5          = false;  ///< SSE5 instructions (AMD)
      bool avx           = false;  ///< AVX instructions
      bool avx2          = false;  ///< AVX2 instructions
      bool avx512f       = false;  ///< AVX-512 Foundation instructions
      bool bmi1          = false;  ///< Bit Manipulation Instruction Set 1
      bool bmi2          = false;  ///< Bit Manipulation Instruction Set 2
      bool fma           = false;  ///< Fused Multiply-Add
      bool aes           = false;  ///< AES instruction set
      bool sha           = false;  ///< SHA extensions
      bool x64           = false;  ///< 64-bit mode (Long Mode)
      bool rdrand        = false;  ///< RDRAND instruction
      bool rdseed        = false;  ///< RDSEED instruction
      bool adx           = false;  ///< Multi-precision Add-Carry
      bool clmul         = false;  ///< Carry-less Multiplication
      bool pclmulqdq     = false;  ///< PCLMULQDQ instruction
      bool smap          = false;  ///< Supervisor Mode Access Prevention
      bool smep          = false;  ///< Supervisor Mode Execution Protection
      bool tsc           = false;  ///< Time Stamp Counter
      bool tscinvariant  = false;  ///< Invariant TSC
      bool invariant_tsc = false;  ///< Invariant Time Stamp Counter
      bool constant_tsc  = false;  ///< Constant Time Stamp Counter
    };

    /**
     * @brief Queries x86/x64 CPU features via CPUID instruction.
     * @return CPUFeatures structure with detected features.
     *
     * @note Only valid for x86/x64 architectures. Returns empty structure
     *       for other architectures.
     */
    static CPUFeatures get_features()
    {
      CPUFeatures features;

      if ( get_architecture() != ArchType::X86 && get_architecture() != ArchType::X64 ) { return features; }

      CPUIDResult cpuid1, cpuid7, cpuidExt;

      // Basic CPUID level 1
      cpuId( cpuid1, 1 );

      // Check standard features
      features.mmx       = ( cpuid1.edx & EDX_MMX_bit ) != 0;
      features.sse       = ( cpuid1.edx & EDX_SSE_bit ) != 0;
      features.sse2      = ( cpuid1.edx & EDX_SSE2_bit ) != 0;
      features.sse3      = ( cpuid1.ecx & ECX_SSE3_bit ) != 0;
      features.ssse3     = ( cpuid1.ecx & ECX_SSSE3_bit ) != 0;
      features.sse41     = ( cpuid1.ecx & ECX_SSE41_bit ) != 0;
      features.sse42     = ( cpuid1.ecx & ECX_SSE42_bit ) != 0;
      features.avx       = ( cpuid1.ecx & ECX_AVX_bit ) != 0;
      features.x64       = ( cpuid1.edx & ( 1UL << 29 ) ) != 0;  // LM bit
      features.tsc       = ( cpuid1.edx & ( 1UL << 4 ) ) != 0;
      features.rdrand    = ( cpuid1.ecx & ( 1UL << 30 ) ) != 0;
      features.aes       = ( cpuid1.ecx & ( 1UL << 25 ) ) != 0;
      features.pclmulqdq = ( cpuid1.ecx & ( 1UL << 1 ) ) != 0;
      features.fma       = ( cpuid1.ecx & ( 1UL << 12 ) ) != 0;

      // Check extended features (level 7)
      cpuId( cpuid7, 7, 0 );
      features.avx2    = ( cpuid7.ebx & ECX_AVX2_bit ) != 0;
      features.avx512f = ( cpuid7.ebx & ECX_AVX512F_bit ) != 0;
      features.bmi1    = ( cpuid7.ebx & ( 1UL << 3 ) ) != 0;
      features.bmi2    = ( cpuid7.ebx & ( 1UL << 8 ) ) != 0;
      features.sha     = ( cpuid7.ebx & ( 1UL << 29 ) ) != 0;
      features.rdseed  = ( cpuid7.ebx & ( 1UL << 18 ) ) != 0;
      features.adx     = ( cpuid7.ebx & ( 1UL << 19 ) ) != 0;
      features.smap    = ( cpuid7.ebx & ( 1UL << 20 ) ) != 0;
      features.smep    = ( cpuid7.ebx & ( 1UL << 7 ) ) != 0;
      features.clmul   = ( cpuid7.ebx & ( 1UL << 1 ) ) != 0;

      // Extended CPUID for AMD-specific features
      CPUIDResult cpuidExtMax;
      cpuId( cpuidExtMax, 0x80000000 );
      if ( cpuidExtMax.eax >= 0x80000001 )
      {
        cpuId( cpuidExt, 0x80000001 );
        features.mmxplus = ( cpuidExt.edx & EDX_MMXplus_bit ) != 0;
        features.sse4a   = ( cpuidExt.ecx & ECX_SSE4A_bit ) != 0;
        features.sse5    = ( cpuidExt.ecx & ECX_SSE5_bit ) != 0;
      }

      // Check TSC features
      if ( cpuid1.edx & ( 1UL << 4 ) )
      {  // TSC
        CPUIDResult cpuid80000007;
        if ( cpuidExtMax.eax >= 0x80000007 )
        {
          cpuId( cpuid80000007, 0x80000007 );
          features.invariant_tsc = ( cpuid80000007.edx & ( 1UL << 8 ) ) != 0;
        }
      }

      return features;
    }

    /// @name x86/x64 Feature Checkers
    /// @brief Convenience functions to check specific CPU features.
    /// @{
    static bool has_mmx() { return get_features().mmx; }
    static bool has_sse() { return get_features().sse; }
    static bool has_sse2() { return get_features().sse2; }
    static bool has_sse3() { return get_features().sse3; }
    static bool has_ssse3() { return get_features().ssse3; }
    static bool has_sse41() { return get_features().sse41; }
    static bool has_sse42() { return get_features().sse42; }
    static bool has_avx() { return get_features().avx; }
    static bool has_avx2() { return get_features().avx2; }
    static bool has_avx512f() { return get_features().avx512f; }
    static bool has_fma() { return get_features().fma; }
    static bool has_aes() { return get_features().aes; }
    /// @}

#endif  // UTILS_CPU_X86

    // =========================================================================
    // ARM CPU DETECTION
    // =========================================================================

#ifdef UTILS_CPU_ARM

    /**
     * @struct ARMFeatures
     * @brief Holds ARM CPU feature flags.
     */
    struct ARMFeatures
    {
      bool neon    = false;  ///< NEON SIMD instructions
      bool crypto  = false;  ///< Cryptographic extensions
      bool crc32   = false;  ///< CRC32 instructions
      bool sha1    = false;  ///< SHA1 instructions
      bool sha2    = false;  ///< SHA2 instructions
      bool aes     = false;  ///< AES instructions
      bool pmull   = false;  ///< Polynomial Multiply Long
      bool fp16    = false;  ///< Half-precision floating point
      bool dotprod = false;  ///< Dot Product instructions
      bool sve     = false;  ///< Scalable Vector Extension
      bool sve2    = false;  ///< Scalable Vector Extension 2
      bool bf16    = false;  ///< Brain Float 16 support
      bool i8mm    = false;  ///< Int8 Matrix Multiply
      bool fhm     = false;  ///< FP16 Fused Multiply-Add
    };

    /**
     * @brief Queries ARM CPU features.
     * @return ARMFeatures structure with detected features.
     *
     * On Linux, reads /proc/cpuinfo. On Windows, uses IsProcessorFeaturePresent.
     */
    static ARMFeatures get_arm_features()
    {
      ARMFeatures features;

#if defined( __linux__ ) && ( defined( __arm__ ) || defined( __aarch64__ ) )
      // Read CPU features from /proc/cpuinfo on Linux
      FILE * cpuinfo = fopen( "/proc/cpuinfo", "r" );
      if ( cpuinfo )
      {
        char line[256];
        while ( fgets( line, sizeof( line ), cpuinfo ) )
        {
          if ( strstr( line, "Features" ) )
          {
            features.neon    = strstr( line, "neon" ) != nullptr;
            features.crypto  = strstr( line, "crypto" ) != nullptr;
            features.crc32   = strstr( line, "crc32" ) != nullptr;
            features.sha1    = strstr( line, "sha1" ) != nullptr;
            features.sha2    = strstr( line, "sha2" ) != nullptr;
            features.aes     = strstr( line, "aes" ) != nullptr;
            features.pmull   = strstr( line, "pmull" ) != nullptr;
            features.fp16    = strstr( line, "fp16" ) != nullptr;
            features.dotprod = strstr( line, "dotprod" ) != nullptr;
            features.sve     = strstr( line, "sve" ) != nullptr;
            features.sve2    = strstr( line, "sve2" ) != nullptr;
            features.bf16    = strstr( line, "bf16" ) != nullptr;
            features.i8mm    = strstr( line, "i8mm" ) != nullptr;
            features.fhm     = strstr( line, "fhm" ) != nullptr;
            break;
          }
        }
        fclose( cpuinfo );
      }
#elif defined( _WIN32 ) && ( defined( _M_ARM ) || defined( _M_ARM64 ) )
      // Windows ARM feature detection
      features.neon = true;  // All Windows on ARM devices support NEON

      // Check for crypto extensions
      typedef BOOL( WINAPI * IsProcessorFeaturePresentFunc )( DWORD );
      HMODULE kernel32 = GetModuleHandleA( "kernel32.dll" );
      if ( kernel32 )
      {
        auto pIsProcessorFeaturePresent = (IsProcessorFeaturePresentFunc) GetProcAddress(
          kernel32,
          "IsProcessorFeaturePresent" );
        if ( pIsProcessorFeaturePresent )
        {
          features.crypto = pIsProcessorFeaturePresent( PF_ARM_V8_CRYPTO_INSTRUCTIONS_AVAILABLE ) != 0;
          features.crc32  = pIsProcessorFeaturePresent( PF_ARM_V8_CRC32_INSTRUCTIONS_AVAILABLE ) != 0;
        }
      }
#endif

      return features;
    }

    /// @name ARM Feature Checkers
    /// @brief Convenience functions to check specific ARM CPU features.
    /// @{
    static bool has_neon() { return get_arm_features().neon; }
    static bool has_arm_crypto() { return get_arm_features().crypto; }
    static bool has_arm_crc32() { return get_arm_features().crc32; }
    static bool has_arm_aes() { return get_arm_features().aes; }
    /// @}

#endif  // UTILS_CPU_ARM

    // =========================================================================
    // POWERPC CPU DETECTION
    // =========================================================================

#ifdef UTILS_CPU_PPC

    /**
     * @struct PowerPCFeatures
     * @brief Holds PowerPC CPU feature flags.
     */
    struct PowerPCFeatures
    {
      bool altivec = false;  ///< AltiVec SIMD instructions
      bool vsx     = false;  ///< Vector Scalar Extension
      bool vsx2    = false;  ///< VSX version 2
      bool vsx3    = false;  ///< VSX version 3
      bool power8  = false;  ///< POWER8 architecture
      bool power9  = false;  ///< POWER9 architecture
      bool power10 = false;  ///< POWER10 architecture
    };

    /**
     * @brief Queries PowerPC CPU features.
     * @return PowerPCFeatures structure with detected features.
     *
     * Reads /proc/cpuinfo on Linux PowerPC systems.
     */
    static PowerPCFeatures get_ppc_features()
    {
      PowerPCFeatures features;

#if defined( __linux__ ) && ( defined( __powerpc__ ) || defined( __ppc__ ) )
      FILE * cpuinfo = fopen( "/proc/cpuinfo", "r" );
      if ( cpuinfo )
      {
        char line[256];
        while ( fgets( line, sizeof( line ), cpuinfo ) )
        {
          if ( strstr( line, "cpu" ) )
          {
            features.altivec = strstr( line, "altivec" ) != nullptr;
            features.vsx     = strstr( line, "vsx" ) != nullptr;

            // Check for Power8+ features
            if ( strstr( line, "POWER8" ) || strstr( line, "POWER9" ) || strstr( line, "POWER10" ) )
            {
              features.power8 = true;
              if ( strstr( line, "POWER9" ) || strstr( line, "POWER10" ) )
              {
                features.power9 = true;
                if ( strstr( line, "POWER10" ) ) { features.power10 = true; }
              }
            }
            break;
          }
        }
        fclose( cpuinfo );
      }
#endif

      return features;
    }

    /// @brief Checks if AltiVec instructions are available.
    static bool has_altivec() { return get_ppc_features().altivec; }
    /// @brief Checks if VSX instructions are available.
    static bool has_vsx() { return get_ppc_features().vsx; }

#endif  // UTILS_CPU_PPC

    // =========================================================================
    // MIPS CPU DETECTION
    // =========================================================================

#ifdef UTILS_CPU_MIPS

    /**
     * @struct MIPSFeatures
     * @brief Holds MIPS CPU feature flags.
     */
    struct MIPSFeatures
    {
      bool msa   = false;  ///< MIPS SIMD Architecture
      bool dspr2 = false;  ///< DSP Revision 2
      bool dspr3 = false;  ///< DSP Revision 3
      bool eva   = false;  ///< Enhanced Virtual Addressing
    };

    /**
     * @brief Queries MIPS CPU features.
     * @return MIPSFeatures structure with detected features.
     *
     * Reads /proc/cpuinfo on Linux MIPS systems.
     */
    static MIPSFeatures get_mips_features()
    {
      MIPSFeatures features;

#if defined( __linux__ ) && ( defined( __mips__ ) || defined( __mips64__ ) )
      FILE * cpuinfo = fopen( "/proc/cpuinfo", "r" );
      if ( cpuinfo )
      {
        char line[256];
        while ( fgets( line, sizeof( line ), cpuinfo ) )
        {
          if ( strstr( line, "ASEs implemented" ) )
          {
            features.msa   = strstr( line, "msa" ) != nullptr;
            features.dspr2 = strstr( line, "dsp2" ) != nullptr;
            features.dspr3 = strstr( line, "dsp3" ) != nullptr;
            features.eva   = strstr( line, "eva" ) != nullptr;
            break;
          }
        }
        fclose( cpuinfo );
      }
#endif

      return features;
    }

    /// @brief Checks if MSA (MIPS SIMD Architecture) is available.
    static bool has_msa() { return get_mips_features().msa; }

#endif  // UTILS_CPU_MIPS

    // =========================================================================
    // RISC-V CPU DETECTION
    // =========================================================================

#ifdef UTILS_CPU_RISCV

    /**
     * @struct RISCVFeatures
     * @brief Holds RISC-V CPU feature flags.
     */
    struct RISCVFeatures
    {
      bool rva20u64    = false;  ///< RV64GC (general purpose 64-bit)
      bool vector      = false;  ///< Vector extension (V)
      bool packed_simd = false;  ///< Packed SIMD extension (P)
      bool bitmanip    = false;  ///< Bit manipulation extension (B)
      bool crypto      = false;  ///< Cryptographic extension (K)
    };

    /**
     * @brief Queries RISC-V CPU features.
     * @return RISCVFeatures structure with detected features.
     *
     * Reads /proc/cpuinfo on Linux RISC-V systems.
     */
    static RISCVFeatures get_riscv_features()
    {
      RISCVFeatures features;

#if defined( __linux__ ) && defined( __riscv )
      FILE * cpuinfo = fopen( "/proc/cpuinfo", "r" );
      if ( cpuinfo )
      {
        char line[256];
        while ( fgets( line, sizeof( line ), cpuinfo ) )
        {
          if ( strstr( line, "isa" ) )
          {
            features.rva20u64    = strstr( line, "rv64gc" ) != nullptr || strstr( line, "rv64imafdc" ) != nullptr;
            features.vector      = strstr( line, "v" ) != nullptr;
            features.packed_simd = strstr( line, "p" ) != nullptr;
            features.bitmanip    = strstr( line, "b" ) != nullptr;
            features.crypto      = strstr( line, "k" ) != nullptr;
            break;
          }
        }
        fclose( cpuinfo );
      }
#endif

      return features;
    }

#endif  // UTILS_CPU_RISCV

    // =========================================================================
    // COMMON CPU INFORMATION
    // =========================================================================

    /**
     * @brief Gets the CPU vendor string.
     * @return Vendor string (e.g., "GenuineIntel", "AuthenticAMD").
     *
     * For x86/x64, uses CPUID. For other architectures, reads from /proc/cpuinfo
     * on Linux systems.
     */
    static std::string get_cpu_vendor()
    {
      switch ( get_architecture() )
      {
#ifdef UTILS_CPU_X86
        case ArchType::X86:
        case ArchType::X64:
        {
          CPUIDResult cpuid0;
          cpuId( cpuid0, 0 );
          char vendor[13] = { 0 };
          memcpy( vendor, &cpuid0.ebx, 4 );
          memcpy( vendor + 4, &cpuid0.edx, 4 );
          memcpy( vendor + 8, &cpuid0.ecx, 4 );
          return vendor;
        }
#endif
        default:
        {
#if defined( __linux__ )
          FILE * cpuinfo = fopen( "/proc/cpuinfo", "r" );
          if ( cpuinfo )
          {
            char line[256];
            while ( fgets( line, sizeof( line ), cpuinfo ) )
            {
              if ( strstr( line, "vendor_id" ) || strstr( line, "CPU implementer" ) )
              {
                fclose( cpuinfo );
                char * colon = strchr( line, ':' );
                if ( colon )
                {
                  colon++;
                  while ( *colon == ' ' || *colon == '\t' ) colon++;
                  // Remove trailing newline
                  char * nl = strchr( colon, '\n' );
                  if ( nl ) *nl = 0;
                  return colon;
                }
              }
            }
            fclose( cpuinfo );
          }
#endif
          return "Unknown";
        }
      }
    }

    /**
     * @brief Gets the CPU model name string.
     * @return Model name string (e.g., "Intel(R) Core(TM) i7-10700K CPU @ 3.80GHz").
     *
     * Reads from /proc/cpuinfo on Linux systems.
     */
    static std::string get_cpu_model()
    {
#if defined( __linux__ )
      FILE * cpuinfo = fopen( "/proc/cpuinfo", "r" );
      if ( cpuinfo )
      {
        char line[256];
        while ( fgets( line, sizeof( line ), cpuinfo ) )
        {
          if ( strstr( line, "model name" ) || strstr( line, "Processor" ) || strstr( line, "cpu model" ) )
          {
            fclose( cpuinfo );
            char * colon = strchr( line, ':' );
            if ( colon )
            {
              colon++;
              while ( *colon == ' ' || *colon == '\t' ) colon++;
              char * nl = strchr( colon, '\n' );
              if ( nl ) *nl = 0;
              return colon;
            }
          }
        }
        fclose( cpuinfo );
      }
#endif
      return "Unknown";
    }

    /**
     * @brief Gets the number of available CPU cores.
     * @return Number of logical CPU cores.
     *
     * Uses system-specific APIs: GetSystemInfo (Windows), sysconf (Linux/macOS).
     */
    static int get_cpu_count()
    {
#if defined( _WIN32 ) || defined( _WIN64 )
      SYSTEM_INFO sysInfo;
      GetSystemInfo( &sysInfo );
      return sysInfo.dwNumberOfProcessors;
#elif defined( __linux__ ) || defined( __APPLE__ )
      return sysconf( _SC_NPROCESSORS_ONLN );
#else
      return 1;
#endif
    }

    /**
     * @brief Gets a formatted string with comprehensive CPU information.
     * @return Formatted string containing architecture, vendor, model, and core count.
     */
    static std::string get_cpu_info_string()
    {
      std::string arch   = get_architecture_string();
      std::string vendor = get_cpu_vendor();
      std::string model  = get_cpu_model();
      int         cores  = get_cpu_count();

      return fmt::format( "Architecture: {}, Vendor: {}, Model: {}, Cores: {}", arch, vendor, model, cores );
    }
  };

  /*==========================================================================*\
  |                            SYSTEM UTILITIES                                |
  \*==========================================================================*/

  /**
   * @brief Extracts the directory portion from a full path.
   * @param path Full file system path.
   * @return Directory portion of the path.
   */
  inline std::string get_basename( std::string_view path )
  {
    namespace fs = std::filesystem;
    return fs::path( path ).parent_path().string();
  }

  /**
   * @brief Extracts the filename portion from a full path.
   * @param path Full file system path.
   * @return Filename portion of the path.
   */
  inline std::string get_filename( std::string_view path )
  {
    namespace fs = std::filesystem;
    return fs::path( path ).filename().string();
  }

  /**
   * @brief Extracts the file extension from a path.
   * @param path Full file system path.
   * @return File extension (including the dot).
   */
  inline std::string get_extension( std::string_view path )
  {
    namespace fs = std::filesystem;
    return fs::path( path ).extension().string();
  }

  /*==========================================================================*\
  |                            OS-SPECIFIC IMPLEMENTATIONS                     |
  \*==========================================================================*/

#if defined( _WIN32 ) || defined( _WIN64 )

#include <windows.h>
#include <winsock2.h>
#include <iphlpapi.h>
#include <direct.h>
#include <ws2tcpip.h>

#pragma comment( lib, "Ws2_32.lib" )
#pragma comment( lib, "Iphlpapi.lib" )

  /**
   * @brief Initializes Windows Sockets (Winsock) API.
   *
   * Uses RAII pattern to ensure cleanup on program exit.
   */
  inline void init_winsock()
  {
    static struct WinSockInit
    {
      WinSockInit()
      {
        WORD    wVersionRequested = MAKEWORD( 2, 2 );
        WSADATA wsaData;
        WSAStartup( wVersionRequested, &wsaData );
      }
      ~WinSockInit() { WSACleanup(); }
    } winsock_init;
    (void) winsock_init;
  }

  /**
   * @brief Gets the value of an environment variable (Windows).
   * @param ename Environment variable name.
   * @param[out] res Value of the environment variable.
   * @return true if variable exists, false otherwise.
   */
  inline bool get_environment( std::string_view ename, std::string & res )
  {
    char  buffer[1024];
    DWORD var_size = GetEnvironmentVariableA( ename.data(), buffer, 1024 );
    if ( var_size == 0 ) return false;
    res = std::string{ buffer, var_size };
    return true;
  }

  /**
   * @brief Sets an environment variable (Windows).
   * @param ename Environment variable name.
   * @param newval New value to set.
   * @param overwrite Unused on Windows (always overwrites).
   */
  inline void set_environment( std::string_view ename, std::string_view newval, bool /* overwrite */ )
  {
    SetEnvironmentVariableA( ename.data(), newval.data() );
  }

  /**
   * @brief Gets the hostname of the system (Windows).
   * @return Hostname string, or empty string on error.
   */
  inline std::string get_host_name()
  {
    char szHostName[1024];
    if ( gethostname( szHostName, 1024 ) == 0 ) { return std::string( szHostName ); }
    return std::string{ "" };
  }

  /**
   * @brief Gets MAC addresses of network interfaces (Windows).
   * @param[out] addr Map of interface names to MAC addresses.
   */
  inline void get_MAC_address( std::map<std::string, std::string> & addr )
  {
    init_winsock();

    IP_ADAPTER_INFO AdapterInfo[16];
    DWORD           dwBufLen = sizeof( AdapterInfo );
    DWORD           dwStatus = GetAdaptersInfo( AdapterInfo, &dwBufLen );

    if ( dwStatus != ERROR_SUCCESS )
    {
      // UTILS_ASSERT0 removed for simplicity
      return;
    }

    PIP_ADAPTER_INFO pAdapterInfo = AdapterInfo;
    addr.clear();

    do
    {
      if ( pAdapterInfo->Type == MIB_IF_TYPE_ETHERNET || pAdapterInfo->Type == IF_TYPE_IEEE80211 )
      {
        unsigned char const * MACData{ pAdapterInfo->Address };
        std::string           str = fmt::format(
          "{:02x}:{:02x}:{:02x}:{:02x}:{:02x}:{:02x}",
          MACData[0],
          MACData[1],
          MACData[2],
          MACData[3],
          MACData[4],
          MACData[5] );
        addr[pAdapterInfo->AdapterName] = str;
      }
      pAdapterInfo = pAdapterInfo->Next;
    } while ( pAdapterInfo );
  }

  /**
   * @brief Gets IP addresses of network interfaces (Windows).
   * @param[out] addr Vector of IP address strings.
   */
  inline void get_IP_address( std::vector<std::string> & addr )
  {
    init_winsock();

    IP_ADAPTER_INFO AdapterInfo[16];
    DWORD           dwBufLen = sizeof( AdapterInfo );
    DWORD           dwStatus = GetAdaptersInfo( AdapterInfo, &dwBufLen );

    if ( dwStatus != ERROR_SUCCESS ) { return; }

    PIP_ADAPTER_INFO pAdapterInfo = AdapterInfo;
    addr.clear();

    do
    {
      addr.emplace_back( pAdapterInfo->IpAddressList.IpAddress.String );
      pAdapterInfo = pAdapterInfo->Next;
    } while ( pAdapterInfo );
  }

  /**
   * @brief Gets the current username (Windows).
   * @return Username string, or empty string on error.
   */
  inline std::string get_user_name()
  {
    char  buffer[1024] = { 0 };
    DWORD size         = 1024;
    if ( GetUserNameA( buffer, &size ) ) { return std::string{ buffer }; }
    return std::string{ "" };
  }

  /**
   * @brief Gets the home directory path (Windows).
   * @return Home directory path, or empty string on error.
   */
  inline std::string get_home_directory()
  {
    char buffer[1024];
    if ( GetEnvironmentVariableA( "USERPROFILE", buffer, 1024 ) ) { return std::string{ buffer }; }
    return std::string{ "" };
  }

  /**
   * @brief Gets the full path of the current executable (Windows).
   * @return Executable path string.
   */
  inline std::string get_executable_path_name()
  {
    char  pathName[1024];
    DWORD pathNameCapacity = 1024;
    GetModuleFileNameA( nullptr, pathName, pathNameCapacity );
    return std::string{ pathName };
  }

  /**
   * @brief Checks if a file exists (Windows).
   * @param fname Path to the file.
   * @return true if file exists and is a regular file, false otherwise.
   */
  inline bool check_if_file_exists( std::string_view fname )
  {
    DWORD ftyp = GetFileAttributesA( fname.data() );
    if ( ftyp == INVALID_FILE_ATTRIBUTES ) return false;
    if ( ftyp & FILE_ATTRIBUTE_DIRECTORY ) return false;
    return true;
  }

  /**
   * @brief Checks if a directory exists (Windows).
   * @param dirname Path to the directory.
   * @return true if directory exists, false otherwise.
   */
  inline bool check_if_dir_exists( std::string_view dirname )
  {
    DWORD ftyp = GetFileAttributesA( dirname.data() );
    if ( ftyp == INVALID_FILE_ATTRIBUTES ) return false;
    if ( ftyp & FILE_ATTRIBUTE_DIRECTORY ) return true;
    return false;
  }

  /**
   * @brief Creates a directory (Windows).
   * @param dirname Path of the directory to create.
   * @param mode Unused on Windows (ignored).
   * @return true if directory created or already exists, false on error.
   */
  inline bool make_directory( std::string_view dirname, unsigned /* mode */ = 0777 )
  {
    if ( check_if_dir_exists( dirname ) ) return true;
    return CreateDirectoryA( dirname.data(), nullptr ) != 0;
  }

  /**
   * @brief Gets the current date in YYYY-MM-DD format (Windows).
   * @return Date string.
   */
  inline std::string get_date()
  {
    SYSTEMTIME st;
    GetLocalTime( &st );
    return fmt::format( "{:04}-{:02}-{:02}", st.wYear, st.wMonth, st.wDay );
  }

  /**
   * @brief Gets the current time in HH:MM:SS format (Windows).
   * @return Time string.
   */
  inline std::string get_day_time()
  {
    SYSTEMTIME st;
    GetLocalTime( &st );
    return fmt::format( "{:02}:{:02}:{:02}", st.wHour, st.wMinute, st.wSecond );
  }

  /**
   * @brief Gets a timestamp string for log files (Windows).
   * @return String in format "date_YYYY-MM-DD_time_HH-MM-SS".
   */
  inline std::string get_log_date_time()
  {
    SYSTEMTIME st;
    GetLocalTime( &st );
    return fmt::format(
      "date_{:04}-{:02}-{:02}_time_{:02}-{:02}-{:02}",
      st.wYear,
      st.wMonth,
      st.wDay,
      st.wHour,
      st.wMinute,
      st.wSecond );
  }

  /**
   * @brief Gets combined date and time (Windows).
   * @return String in format "HH:MM:SS YYYY-MM-DD".
   */
  inline std::string get_day_time_and_date()
  {
    SYSTEMTIME st;
    GetLocalTime( &st );
    return fmt::format(
      "{:02}:{:02}:{:02} {:04}-{:02}-{:02}",
      st.wHour,
      st.wMinute,
      st.wSecond,
      st.wYear,
      st.wMonth,
      st.wDay );
  }

// =========================================================================
// LINUX IMPLEMENTATIONS
// =========================================================================
#elif defined( __linux__ )

#include <arpa/inet.h>
#include <dirent.h>
#include <ifaddrs.h>
#include <net/if.h>
#include <netdb.h>
#include <sys/ioctl.h>
#include <sys/socket.h>
#include <sys/stat.h>
#include <sys/statvfs.h>  // Added for statvfs
#include <unistd.h>
#include <cstdio>
#include <cstring>

  /**
   * @brief Gets the value of an environment variable (Linux).
   * @param ename Environment variable name.
   * @param[out] res Value of the environment variable.
   * @return true if variable exists, false otherwise.
   */
  inline bool get_environment( std::string_view ename, std::string & res )
  {
    char const * RES = getenv( ename.data() );
    if ( RES == nullptr ) return false;
    res = std::string{ RES };
    return true;
  }

  /**
   * @brief Sets an environment variable (Linux).
   * @param ename Environment variable name.
   * @param newval New value to set.
   * @param overwrite If true, overwrite existing variable.
   */
  inline void set_environment( std::string_view ename, std::string_view newval, bool overwrite )
  {
    int res = setenv( ename.data(), newval.data(), overwrite ? 1 : 0 );
    (void) res;  // Suppress unused variable warning
  }

  /**
   * @brief Gets the hostname of the system (Linux).
   * @return Hostname string, or empty string on error.
   */
  inline std::string get_host_name()
  {
    char szHostName[1024];
    bool ok = gethostname( szHostName, 1024 ) == 0;
    if ( ok ) return std::string( szHostName );
    return std::string{ "" };
  }

  /**
   * @brief Gets MAC addresses of Ethernet and wireless interfaces (Linux).
   * @param[out] mac_addr Map of interface names to MAC addresses.
   */
  inline void get_MAC_address( std::map<std::string, std::string> & mac_addr )
  {
    char           buf[8192] = { 0 };
    struct ifconf  ifc       = { 0 };
    struct ifreq * ifr       = nullptr;

    int sck = socket( PF_INET, SOCK_DGRAM, 0 );
    if ( sck < 0 ) return;

    ifc.ifc_len = sizeof( buf );
    ifc.ifc_buf = buf;
    if ( ioctl( sck, SIOCGIFCONF, &ifc ) < 0 )
    {
      close( sck );
      return;
    }

    ifr             = ifc.ifc_req;
    int nInterfaces = ifc.ifc_len / sizeof( struct ifreq );
    if ( nInterfaces <= 0 )
    {
      close( sck );
      return;
    }

    for ( int i = 0; i < nInterfaces; ++i )
    {
      struct ifreq * item = &ifr[i];
      if ( strncmp( "eth", item->ifr_name, 3 ) == 0 || strncmp( "en", item->ifr_name, 2 ) == 0 )
      {
        if ( ioctl( sck, SIOCGIFHWADDR, item ) >= 0 )
        {
          std::string mac_str = fmt::format(
            "{:02x}:{:02x}:{:02x}:{:02x}:{:02x}:{:02x}",
            static_cast<unsigned char>( item->ifr_hwaddr.sa_data[0] ),
            static_cast<unsigned char>( item->ifr_hwaddr.sa_data[1] ),
            static_cast<unsigned char>( item->ifr_hwaddr.sa_data[2] ),
            static_cast<unsigned char>( item->ifr_hwaddr.sa_data[3] ),
            static_cast<unsigned char>( item->ifr_hwaddr.sa_data[4] ),
            static_cast<unsigned char>( item->ifr_hwaddr.sa_data[5] ) );
          mac_addr[item->ifr_name] = mac_str;
        }
      }
    }
    close( sck );
  }

  /**
   * @brief Gets IP addresses of the system (Linux).
   * @param[out] addr Vector of IP address strings.
   */
  inline void get_IP_address( std::vector<std::string> & addr )
  {
    addr.clear();
    char szHostName[128];
    if ( gethostname( szHostName, 128 ) == 0 )
    {
      struct hostent * pHost;
      pHost = gethostbyname( szHostName );

      for ( int i = 0; pHost != nullptr && pHost->h_addr_list[i] != nullptr; i++ )
      {
        unsigned char * h   = reinterpret_cast<unsigned char *>( pHost->h_addr_list[i] );
        std::string     str = "";
        for ( int j = 0; j < pHost->h_length; j++ )
        {
          if ( j > 0 ) str += '.';
          str += fmt::format( "{}", static_cast<int>( h[j] ) );
        }
        addr.emplace_back( str );
      }
    }
  }

  /**
   * @brief Gets the current username (Linux).
   * @return Username string, or empty string if not found.
   */
  inline std::string get_user_name()
  {
    char const * USER = getenv( "USER" );
    if ( USER == nullptr ) return "";
    return std::string{ USER };
  }

  /**
   * @brief Gets the home directory path (Linux).
   * @return Home directory path, or empty string if not found.
   */
  inline std::string get_home_directory()
  {
    char const * HOME = getenv( "HOME" );
    if ( HOME == nullptr ) return "";
    return std::string{ HOME };
  }

  /**
   * @brief Gets the full path of the current executable (Linux).
   * @return Executable path string.
   */
  inline std::string get_executable_path_name()
  {
    char     pathName[1024];
    uint32_t pathNameCapacity = 1024;
    uint32_t pathNameSize     = uint32_t( readlink( "/proc/self/exe", pathName, pathNameCapacity - 1 ) );
    if ( pathNameSize > 0 )
    {
      pathName[pathNameSize] = '\0';
      return std::string{ pathName };
    }
    return std::string{ "" };
  }

  /**
   * @brief Checks if a file exists (Linux).
   * @param fname Path to the file.
   * @return true if file exists and is a regular file, false otherwise.
   */
  inline bool check_if_file_exists( std::string_view fname )
  {
    struct stat buffer;
    if ( stat( fname.data(), &buffer ) == 0 ) return S_ISREG( buffer.st_mode );
    return false;
  }

  /**
   * @brief Checks if a directory exists (Linux).
   * @param dirname Path to the directory.
   * @return true if directory exists, false otherwise.
   */
  inline bool check_if_dir_exists( std::string_view dirname )
  {
    struct stat buffer;
    if ( stat( dirname.data(), &buffer ) == 0 ) return S_ISDIR( buffer.st_mode );
    return false;
  }

  /**
   * @brief Creates a directory (Linux).
   * @param dirname Path of the directory to create.
   * @param mode Permissions for the new directory (octal).
   * @return true if directory created or already exists, false on error.
   */
  inline bool make_directory( std::string_view dirname, unsigned mode = 0777 )
  {
    if ( check_if_dir_exists( dirname ) ) return true;
    return mkdir( dirname.data(), mode_t( mode ) ) == 0;
  }

// =========================================================================
// macOS IMPLEMENTATIONS
// =========================================================================
#elif defined( __APPLE__ )

#include <dirent.h>
#include <ifaddrs.h>
#include <mach-o/dyld.h>
#include <net/if.h>
#include <net/if_dl.h>
#include <netdb.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <cstdio>
#include <cstring>
#include <vector>

  /**
   * @brief Gets the value of an environment variable (macOS).
   * @param ename Environment variable name.
   * @param[out] res Value of the environment variable.
   * @return true if variable exists, false otherwise.
   */
  inline bool get_environment( std::string_view ename, std::string & res )
  {
    char const * RES = getenv( ename.data() );
    if ( RES == nullptr ) return false;
    res = std::string{ RES };
    return true;
  }

  /**
   * @brief Sets an environment variable (macOS).
   * @param ename Environment variable name.
   * @param newval New value to set.
   * @param overwrite If true, overwrite existing variable.
   */
  inline void set_environment( std::string_view ename, std::string_view newval, bool overwrite )
  {
    int res = setenv( ename.data(), newval.data(), overwrite ? 1 : 0 );
    (void) res;  // Suppress unused variable warning
  }

  /**
   * @brief Gets the hostname of the system (macOS).
   * @return Hostname string, or empty string on error.
   */
  inline std::string get_host_name()
  {
    char szHostName[1024];
    bool ok = gethostname( szHostName, 1024 ) == 0;
    if ( ok ) return std::string{ szHostName };
    return std::string{ "" };
  }

  /**
   * @brief Gets MAC addresses of network interfaces (macOS).
   * @param[out] addr Map of interface names to MAC addresses.
   */
  inline void get_MAC_address( std::map<std::string, std::string> & addr )
  {
    struct ifaddrs *ifap, *ifaptr;

    int ier = getifaddrs( &ifap );
    if ( ier < 0 ) return;

    addr.clear();
    for ( ifaptr = ifap; ifaptr != nullptr; ifaptr = ifaptr->ifa_next )
    {
      if ( strncmp( "en", ifaptr->ifa_name, 2 ) == 0 && ifaptr->ifa_addr->sa_family == AF_LINK )
      {
        struct sockaddr_dl * tmp = reinterpret_cast<struct sockaddr_dl *>( ( ifaptr )->ifa_addr );
        unsigned char *      ptr = reinterpret_cast<unsigned char *>( LLADDR( tmp ) );

        std::string mac_str =
          fmt::format( "{:02x}:{:02x}:{:02x}:{:02x}:{:02x}:{:02x}", ptr[0], ptr[1], ptr[2], ptr[3], ptr[4], ptr[5] );
        addr[ifaptr->ifa_name] = mac_str;
      }
    }
    freeifaddrs( ifap );
  }

  /**
   * @brief Gets IP addresses of the system (macOS).
   * @param[out] addr Vector of IP address strings.
   */
  inline void get_IP_address( std::vector<std::string> & addr )
  {
    addr.clear();
    char szHostName[128];
    if ( gethostname( szHostName, 128 ) == 0 )
    {
      struct hostent * pHost;
      pHost = gethostbyname( szHostName );

      for ( int i = 0; pHost != nullptr && pHost->h_addr_list[i] != nullptr; i++ )
      {
        unsigned char * h = reinterpret_cast<unsigned char *>( pHost->h_addr_list[i] );
        std::string     str{ "" };
        for ( int j = 0; j < pHost->h_length; j++ )
        {
          if ( j > 0 ) str += '.';
          str += fmt::format( "{}", static_cast<int>( h[j] ) );
        }
        addr.emplace_back( str );
      }
    }
  }

  /**
   * @brief Gets the current username (macOS).
   * @return Username string, or empty string if not found.
   */
  inline std::string get_user_name()
  {
    char const * USER = getenv( "USER" );
    if ( USER == nullptr ) return "";
    return std::string{ USER };
  }

  /**
   * @brief Gets the home directory path (macOS).
   * @return Home directory path, or empty string if not found.
   */
  inline std::string get_home_directory()
  {
    char const * HOME = getenv( "HOME" );
    if ( HOME == nullptr ) return "";
    return std::string{ HOME };
  }

  /**
   * @brief Gets the full path of the current executable (macOS).
   * @return Executable path string.
   */
  inline std::string get_executable_path_name()
  {
    uint32_t pathNameSize = 0;
    _NSGetExecutablePath( nullptr, &pathNameSize );
    std::vector<char> res( pathNameSize + 1 );
    std::fill( res.begin(), res.end(), '\0' );
    if ( _NSGetExecutablePath( res.data(), &pathNameSize ) != 0 ) return "NOT FOUND";
    return std::string{ res.data() };
  }

  /**
   * @brief Checks if a file exists (macOS).
   * @param fname Path to the file.
   * @return true if file exists, false otherwise.
   */
  inline bool check_if_file_exists( std::string_view fname )
  {
    struct stat buffer;
    return ( stat( fname.data(), &buffer ) == 0 );
  }

  /**
   * @brief Checks if a directory exists (macOS).
   * @param dirname Path to the directory.
   * @return true if directory exists, false otherwise.
   */
  inline bool check_if_dir_exists( std::string_view dirname )
  {
    DIR * pDir    = opendir( dirname.data() );
    bool  bExists = false;
    if ( pDir != nullptr )
    {
      bExists = true;
      (void) closedir( pDir );
    }
    return bExists;
  }

  /**
   * @brief Creates a directory (macOS).
   * @param dirname Path of the directory to create.
   * @param mode Permissions for the new directory (octal).
   * @return true if directory created or already exists, false on error.
   */
  inline bool make_directory( std::string_view dirname, unsigned mode = 0777 )
  {
    if ( check_if_dir_exists( dirname ) ) return true;
    return mkdir( dirname.data(), mode_t( mode ) ) == 0;
  }

#else
#error "Unsupported operating system!"
#endif

// =========================================================================
// COMMON FUNCTIONS FOR LINUX AND macOS (date and time)
// =========================================================================
#if defined( __linux__ ) || defined( __APPLE__ )

  /**
   * @brief Gets the current date in YYYY-MM-DD format (Linux/macOS).
   * @return Date string.
   */
  inline std::string get_date()
  {
    char   buffer[20];
    time_t rawtime;
    time( &rawtime );
    tm timeinfo;
    localtime_r( &rawtime, &timeinfo );
    strftime( buffer, 20, "%F", &timeinfo );
    return std::string{ buffer };
  }

  /**
   * @brief Gets the current time in HH:MM:SS format (Linux/macOS).
   * @return Time string.
   */
  inline std::string get_day_time()
  {
    char   buffer[20];
    time_t rawtime;
    time( &rawtime );
    tm timeinfo;
    localtime_r( &rawtime, &timeinfo );
    strftime( buffer, 20, "%T", &timeinfo );
    return std::string{ buffer };
  }

  /**
   * @brief Gets combined date and time (Linux/macOS).
   * @return String in format "HH:MM:SS YYYY-MM-DD".
   */
  inline std::string get_day_time_and_date()
  {
    char   buffer[40];
    time_t rawtime;
    time( &rawtime );
    tm timeinfo;
    localtime_r( &rawtime, &timeinfo );
    strftime( buffer, 40, "%T %F", &timeinfo );
    return std::string{ buffer };
  }

  /**
   * @brief Gets a timestamp string for log files (Linux/macOS).
   * @return String in format "date_YYYY-MM-DD_time_HH-MM-SS".
   */
  inline std::string get_log_date_time()
  {
    char   buffer[100];
    time_t rawtime;
    time( &rawtime );
    tm timeinfo;
    localtime_r( &rawtime, &timeinfo );
    strftime( buffer, 100, "date_%Y-%m-%d_time_%H-%M-%S", &timeinfo );
    return std::string{ buffer };
  }

#endif

  /*==========================================================================*\
  |                            SYSTEM INFORMATION                              |
  \*==========================================================================*/

  /**
   * @struct DiskSpaceInfo
   * @brief Holds disk space statistics for a filesystem.
   */
  struct DiskSpaceInfo
  {
    uint64_t total_bytes      = 0;    ///< Total disk space in bytes
    uint64_t free_bytes       = 0;    ///< Free disk space in bytes
    uint64_t used_bytes       = 0;    ///< Used disk space in bytes
    double   usage_percentage = 0.0;  ///< Disk usage percentage (0-100)

    /**
     * @brief Returns a formatted string representation.
     * @return String with human-readable disk space information.
     */
    std::string to_string() const
    {
      return fmt::format(
        "Total: {:.2f} GB, Free: {:.2f} GB, Used: {:.2f} GB ({:.1f}%)",
        total_bytes / ( 1024.0 * 1024.0 * 1024.0 ),
        free_bytes / ( 1024.0 * 1024.0 * 1024.0 ),
        used_bytes / ( 1024.0 * 1024.0 * 1024.0 ),
        usage_percentage );
    }
  };

  /**
   * @struct MemoryInfo
   * @brief Holds system memory statistics.
   */
  struct MemoryInfo
  {
    uint64_t total_physical     = 0;    ///< Total physical memory in bytes
    uint64_t available_physical = 0;    ///< Available physical memory in bytes
    uint64_t used_physical      = 0;    ///< Used physical memory in bytes
    double   usage_percentage   = 0.0;  ///< Physical memory usage percentage (0-100)

    uint64_t total_virtual     = 0;  ///< Total virtual memory in bytes
    uint64_t available_virtual = 0;  ///< Available virtual memory in bytes
    uint64_t used_virtual      = 0;  ///< Used virtual memory in bytes

    /**
     * @brief Returns a formatted string representation.
     * @return String with human-readable memory information.
     */
    std::string to_string() const
    {
      return fmt::format(
        "Physical: {:.2f} GB total, {:.2f} GB available ({:.1f}% used)",
        total_physical / ( 1024.0 * 1024.0 * 1024.0 ),
        available_physical / ( 1024.0 * 1024.0 * 1024.0 ),
        usage_percentage );
    }
  };

  /**
   * @struct SystemLoadInfo
   * @brief Holds system load average information.
   */
  struct SystemLoadInfo
  {
    double load_1min  = 0.0;  ///< 1-minute load average
    double load_5min  = 0.0;  ///< 5-minute load average
    double load_15min = 0.0;  ///< 15-minute load average

    /**
     * @brief Returns a formatted string representation.
     * @return String with load average information.
     */
    std::string to_string() const
    {
      return fmt::format(
        "Load averages: {:.2f} (1min), {:.2f} (5min), {:.2f} (15min)",
        load_1min,
        load_5min,
        load_15min );
    }
  };

  /**
   * @brief Gets system memory information.
   * @return MemoryInfo structure with current memory statistics.
   *
   * Uses platform-specific APIs: GlobalMemoryStatusEx (Windows),
   * sysinfo (Linux), mach APIs (macOS).
   */
  inline MemoryInfo get_memory_info()
  {
    MemoryInfo info;

#if defined( _WIN32 ) || defined( _WIN64 )
    // Windows memory information
    MEMORYSTATUSEX memStatus;
    memStatus.dwLength = sizeof( memStatus );
    if ( GlobalMemoryStatusEx( &memStatus ) )
    {
      info.total_physical     = memStatus.ullTotalPhys;
      info.available_physical = memStatus.ullAvailPhys;
      info.used_physical      = info.total_physical - info.available_physical;
      info.usage_percentage   = memStatus.dwMemoryLoad;

      info.total_virtual     = memStatus.ullTotalPageFile;
      info.available_virtual = memStatus.ullAvailPageFile;
      info.used_virtual      = info.total_virtual - info.available_virtual;
    }

#elif defined( __linux__ )
    // Linux memory information
    struct sysinfo sysInfo;
    if ( sysinfo( &sysInfo ) == 0 )
    {
      info.total_physical     = sysInfo.totalram * sysInfo.mem_unit;
      info.available_physical = sysInfo.freeram * sysInfo.mem_unit;
      info.used_physical      = info.total_physical - info.available_physical;
      info.usage_percentage   = ( info.used_physical * 100.0 ) / info.total_physical;

      info.total_virtual     = ( sysInfo.totalram + sysInfo.totalswap ) * sysInfo.mem_unit;
      info.available_virtual = ( sysInfo.freeram + sysInfo.freeswap ) * sysInfo.mem_unit;
      info.used_virtual      = info.total_virtual - info.available_virtual;
    }

    // Alternative: read from /proc/meminfo for more detailed information
    FILE * meminfo = fopen( "/proc/meminfo", "r" );
    if ( meminfo )
    {
      char     line[256];
      uint64_t mem_total = 0, mem_free = 0, buffers = 0, cached = 0;

      while ( fgets( line, sizeof( line ), meminfo ) )
      {
        if ( strstr( line, "MemTotal:" ) ) { sscanf( line, "MemTotal: %lu kB", &mem_total ); }
        else if ( strstr( line, "MemFree:" ) ) { sscanf( line, "MemFree: %lu kB", &mem_free ); }
        else if ( strstr( line, "Buffers:" ) ) { sscanf( line, "Buffers: %lu kB", &buffers ); }
        else if ( strstr( line, "Cached:" ) ) { sscanf( line, "Cached: %lu kB", &cached ); }
      }
      fclose( meminfo );

      if ( mem_total > 0 )
      {
        info.total_physical = mem_total * 1024;  // Convert kB to bytes
        // More accurate available memory calculation for Linux
        info.available_physical = ( mem_free + buffers + cached ) * 1024;
        info.used_physical      = info.total_physical - info.available_physical;
        info.usage_percentage   = ( info.used_physical * 100.0 ) / info.total_physical;
      }
    }

#elif defined( __APPLE__ )
    // macOS memory information - using sysctl and mach APIs
    uint64_t mem_size;
    size_t   len = sizeof( mem_size );

    // Get total physical memory using sysctl
    if ( sysctlbyname( "hw.memsize", &mem_size, &len, NULL, 0 ) == 0 )
    {
      info.total_physical = mem_size;

      // Get memory statistics using mach APIs
      vm_size_t              page_size;
      mach_port_t            mach_port = mach_host_self();
      vm_statistics64_data_t vm_stats;
      mach_msg_type_number_t count = HOST_VM_INFO64_COUNT;

      if (
        host_page_size( mach_port, &page_size ) == KERN_SUCCESS &&
        host_statistics64( mach_port, HOST_VM_INFO64, (host_info64_t) &vm_stats, &count ) == KERN_SUCCESS )
      {
        // Calculate available memory (free + inactive)
        uint64_t free_memory     = vm_stats.free_count * page_size;
        uint64_t inactive_memory = vm_stats.inactive_count * page_size;

        info.available_physical = free_memory + inactive_memory;
        info.used_physical      = info.total_physical - info.available_physical;

        if ( info.total_physical > 0 ) { info.usage_percentage = ( info.used_physical * 100.0 ) / info.total_physical; }
      }
    }
#endif

    return info;
  }

  /**
   * @brief Gets disk space information for a specific path.
   * @param path Filesystem path to check.
   * @return DiskSpaceInfo structure with disk statistics for the path.
   *
   * Uses platform-specific APIs: GetDiskFreeSpaceEx (Windows),
   * statvfs (Linux/macOS).
   */
  inline DiskSpaceInfo get_disk_space( const std::string & path )
  {
    DiskSpaceInfo info;

#if defined( _WIN32 ) || defined( _WIN64 )
    // Windows disk space
    ULARGE_INTEGER freeBytesAvailable, totalNumberOfBytes, totalNumberOfFreeBytes;

    if ( GetDiskFreeSpaceExA( path.c_str(), &freeBytesAvailable, &totalNumberOfBytes, &totalNumberOfFreeBytes ) )
    {
      info.total_bytes = totalNumberOfBytes.QuadPart;
      info.free_bytes  = totalNumberOfFreeBytes.QuadPart;
      info.used_bytes  = info.total_bytes - info.free_bytes;

      if ( info.total_bytes > 0 ) { info.usage_percentage = ( info.used_bytes * 100.0 ) / info.total_bytes; }
    }

#elif defined( __linux__ )
    // Linux disk space
    struct statvfs stat;

    if ( statvfs( path.c_str(), &stat ) == 0 )
    {
      info.total_bytes = stat.f_blocks * stat.f_frsize;
      info.free_bytes  = stat.f_bfree * stat.f_frsize;
      info.used_bytes  = info.total_bytes - info.free_bytes;

      if ( info.total_bytes > 0 ) { info.usage_percentage = ( info.used_bytes * 100.0 ) / info.total_bytes; }
    }

#elif defined( __APPLE__ )
    // macOS disk space
    struct statvfs stat;

    if ( statvfs( path.c_str(), &stat ) == 0 )
    {
      info.total_bytes = stat.f_blocks * stat.f_frsize;
      info.free_bytes  = stat.f_bfree * stat.f_frsize;
      info.used_bytes  = info.total_bytes - info.free_bytes;

      if ( info.total_bytes > 0 ) { info.usage_percentage = ( info.used_bytes * 100.0 ) / info.total_bytes; }
    }
#endif

    return info;
  }

  /**
   * @brief Gets system load average information.
   * @return SystemLoadInfo structure with load averages.
   *
   * Uses getloadavg (Linux/macOS) or performance counters (Windows).
   */
  inline SystemLoadInfo get_system_load()
  {
    SystemLoadInfo info;

#if defined( __linux__ )
    // Linux load averages
    double loadavg[3];
    if ( getloadavg( loadavg, 3 ) > 0 )
    {
      info.load_1min  = loadavg[0];
      info.load_5min  = loadavg[1];
      info.load_15min = loadavg[2];
    }

#elif defined( __APPLE__ )
    // macOS load averages - simpler approach
    double loadavg[3];
    if ( getloadavg( loadavg, 3 ) > 0 )
    {
      info.load_1min  = loadavg[0];
      info.load_5min  = loadavg[1];
      info.load_15min = loadavg[2];
    }

#elif defined( _WIN32 ) || defined( _WIN64 )
    // Windows doesn't have load averages in the same way
    // We can use performance counters as an approximation
    PDH_HQUERY           query;
    PDH_HCOUNTER         counter;
    PDH_FMT_COUNTERVALUE value;

    if ( PdhOpenQuery( NULL, 0, &query ) == ERROR_SUCCESS )
    {
      if ( PdhAddCounterA( query, "\\Processor(_Total)\\% Processor Time", 0, &counter ) == ERROR_SUCCESS )
      {
        if ( PdhCollectQueryData( query ) == ERROR_SUCCESS )
        {
          Sleep( 1000 );  // Wait 1 second
          if ( PdhCollectQueryData( query ) == ERROR_SUCCESS )
          {
            if ( PdhGetFormattedCounterValue( counter, PDH_FMT_DOUBLE, NULL, &value ) == ERROR_SUCCESS )
            {
              // Convert CPU usage to a load-like value
              // This is an approximation
              info.load_1min = value.doubleValue / 100.0 * Architecture::get_cpu_count();
            }
          }
        }
        PdhRemoveCounter( counter );
      }
      PdhCloseQuery( query );
    }

    info.load_5min  = info.load_1min;
    info.load_15min = info.load_1min;
#endif

    return info;
  }

  /**
   * @brief Gets disk space information for all mounted filesystems.
   * @return Vector of pairs (mount point, DiskSpaceInfo).
   *
   * Platform-specific implementations:
   * - Windows: Iterates logical drives (A-Z)
   * - Linux: Parses /proc/mounts
   * - macOS: Checks common mount points
   */
  inline std::vector<std::pair<std::string, DiskSpaceInfo>> get_all_disk_spaces()
  {
    std::vector<std::pair<std::string, DiskSpaceInfo>> disks;

#if defined( _WIN32 ) || defined( _WIN64 )
    // Windows drives
    DWORD drives  = GetLogicalDrives();
    char  drive[] = "A:\\";

    for ( int i = 0; i < 26; i++ )
    {
      if ( drives & ( 1 << i ) )
      {
        drive[0] = 'A' + i;

        // Check drive type
        UINT type = GetDriveTypeA( drive );
        if ( type == DRIVE_FIXED || type == DRIVE_REMOVABLE || type == DRIVE_REMOTE )
        {
          DiskSpaceInfo info = get_disk_space( drive );
          disks.emplace_back( drive, info );
        }
      }
    }

#elif defined( __linux__ )
    // Linux mounts from /proc/mounts
    FILE * mounts = fopen( "/proc/mounts", "r" );
    if ( mounts )
    {
      char                  line[512];
      std::set<std::string> processed;  // Avoid duplicates

      while ( fgets( line, sizeof( line ), mounts ) )
      {
        char device[256], mountpoint[256], fstype[256];
        if ( sscanf( line, "%255s %255s %255s", device, mountpoint, fstype ) == 3 )
        {
          // Skip virtual filesystems
          if (
            strcmp( fstype, "tmpfs" ) == 0 || strcmp( fstype, "proc" ) == 0 || strcmp( fstype, "sysfs" ) == 0 ||
            strcmp( fstype, "devtmpfs" ) == 0 || strcmp( fstype, "devpts" ) == 0 || strcmp( fstype, "cgroup" ) == 0 ||
            strcmp( fstype, "cgroup2" ) == 0 || strcmp( fstype, "pstore" ) == 0 || strcmp( fstype, "bpf" ) == 0 ||
            strcmp( fstype, "securityfs" ) == 0 || strcmp( fstype, "debugfs" ) == 0 ||
            strcmp( fstype, "tracefs" ) == 0 || strcmp( fstype, "fusectl" ) == 0 || strcmp( fstype, "configfs" ) == 0 ||
            strcmp( fstype, "hugetlbfs" ) == 0 || strcmp( fstype, "mqueue" ) == 0 || strcmp( fstype, "ramfs" ) == 0 )
          {
            continue;
          }

          // Check if we've already processed this mountpoint
          std::string mp( mountpoint );
          if ( processed.find( mp ) == processed.end() )
          {
            DiskSpaceInfo info = get_disk_space( mp );
            if ( info.total_bytes > 0 )
            {  // Valid filesystem
              disks.emplace_back( mp, info );
              processed.insert( mp );
            }
          }
        }
      }
      fclose( mounts );
    }

#elif defined( __APPLE__ )
    // macOS - simplified approach: just check common mount points
    std::vector<std::string> common_paths = { "/", "/System/Volumes/Data", "/Users", "/Volumes" };

    for ( const auto & path : common_paths )
    {
      DiskSpaceInfo info = get_disk_space( path );
      if ( info.total_bytes > 0 ) { disks.emplace_back( path, info ); }
    }
#endif

    return disks;
  }

  /**
   * @brief Gets system uptime in seconds.
   * @return Uptime in seconds.
   *
   * Platform-specific implementations:
   * - Windows: GetTickCount64
   * - Linux: /proc/uptime
   * - macOS: sysctl KERN_BOOTTIME
   */
  inline uint64_t get_system_uptime()
  {
    uint64_t uptime = 0;

#if defined( _WIN32 ) || defined( _WIN64 )
    // Windows uptime
    uptime = GetTickCount64() / 1000;  // Convert milliseconds to seconds

#elif defined( __linux__ )
    // Linux uptime from /proc/uptime
    FILE * uptime_file = fopen( "/proc/uptime", "r" );
    if ( uptime_file )
    {
      double uptime_seconds;
      if ( fscanf( uptime_file, "%lf", &uptime_seconds ) == 1 ) { uptime = static_cast<uint64_t>( uptime_seconds ); }
      fclose( uptime_file );
    }

#elif defined( __APPLE__ )
    // macOS uptime - using sysctl
    struct timeval boottime;
    size_t         len    = sizeof( boottime );
    int            mib[2] = { CTL_KERN, KERN_BOOTTIME };

    if ( sysctl( mib, 2, &boottime, &len, NULL, 0 ) == 0 )
    {
      time_t now;
      time( &now );
      uptime = static_cast<uint64_t>( now - boottime.tv_sec );
    }
#endif

    return uptime;
  }

  /**
   * @brief Formats uptime seconds into a human-readable string.
   * @param seconds Uptime in seconds.
   * @return Formatted string (e.g., "2 days, 03:45:12" or "03:45:12").
   */
  inline std::string format_uptime( uint64_t seconds )
  {
    uint64_t days = seconds / ( 24 * 3600 );
    seconds %= ( 24 * 3600 );
    uint64_t hours = seconds / 3600;
    seconds %= 3600;
    uint64_t minutes = seconds / 60;
    seconds %= 60;

    if ( days > 0 ) { return fmt::format( "{} days, {:02}:{:02}:{:02}", days, hours, minutes, seconds ); }
    else
    {
      return fmt::format( "{:02}:{:02}:{:02}", hours, minutes, seconds );
    }
  }

  /**
   * @brief Generates a comprehensive system summary string.
   * @return Formatted string with architecture, CPU, memory, load, uptime,
   *         disk space, host information, and current date/time.
   */
  inline std::string get_system_summary()
  {
    fmt::memory_buffer buf;
    auto               out = fmt::appender( buf );

    // Architecture info
    fmt::format_to(
      out,
      "Architecture: {}\n"
      "CPU: {} ({} cores)\n"
      "CPU Vendor: {}\n",
      Architecture::get_architecture_string(),
      Architecture::get_cpu_model(),
      Architecture::get_cpu_count(),
      Architecture::get_cpu_vendor() );

    // Memory info
    MemoryInfo mem = get_memory_info();
    fmt::format_to( out, "Memory: {}\n", mem.to_string() );

    // System load
    SystemLoadInfo load = get_system_load();
    fmt::format_to( out, "System Load: {}\n", load.to_string() );

    // Uptime
    uint64_t uptime = get_system_uptime();
    fmt::format_to( out, "Uptime: {}\n", format_uptime( uptime ) );

    // Disk space (root filesystem)
#if defined( _WIN32 ) || defined( _WIN64 )
    constexpr const char * root_path = "C:\\";
#else
    constexpr const char * root_path = "/";
#endif

    DiskSpaceInfo disk = get_disk_space( root_path );
    fmt::format_to( out, "Root Filesystem: {}\n", disk.to_string() );

    // Host information
    fmt::format_to(
      out,
      "Hostname: {}\n"
      "Username: {}\n"
      "Home Directory: {}\n",
      get_host_name(),
      get_user_name(),
      get_home_directory() );

    // Date and time
    fmt::format_to( out, "Current Date/Time: {}\n", get_day_time_and_date() );

    return fmt::to_string( buf );
  }


}  // namespace Utils

#endif  // SYSTEM_UTILS_HXX
