/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2026                                                      |
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

#ifndef UTILS_TRIDIAGONAL_SOLVER_dot_HH
#define UTILS_TRIDIAGONAL_SOLVER_dot_HH

#include "Utils.hh"
#include "Utils_eigen.hh"
#include <vector>
#include <cassert>
#include <memory>

namespace Utils
{

  // ===========================================================================
  // SOLUTORE TRIDIAGONALE SCALARE
  // ===========================================================================
  template <typename Scalar> class TridiagonalSolver
  {
  public:
    using VecS    = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using MatB    = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using integer = Eigen::Index;

  private:
    integer m_size = 0;
    VecS    m_cprime;
    VecS    m_dprime;

  public:
    TridiagonalSolver() noexcept {}

    explicit TridiagonalSolver( integer n_ ) { resize( n_ ); }

    void resize( integer n_ )
    {
      m_size = n_;
      if ( m_size > 0 )
      {
        if ( m_size > 1 ) m_cprime.resize( m_size - 1 );
        m_dprime.resize( m_size );
      }
    }

    integer size() const noexcept { return m_size; }

    /*============================================================
      FACTORIZATION - OTTIMIZZATA
    ============================================================*/

    // Versione template unica - nessun overload
    template <typename DerivedA, typename DerivedB, typename DerivedC> void factorize(
      const Eigen::DenseBase<DerivedA> & a,
      const Eigen::DenseBase<DerivedB> & b,
      const Eigen::DenseBase<DerivedC> & c )
    {
      assert(
        static_cast<integer>( a.size() ) == m_size - 1 && static_cast<integer>( b.size() ) == m_size &&
        static_cast<integer>( c.size() ) == m_size - 1 );

      if ( m_size < 2 ) return;

      // Accesso diretto ai dati per migliori performance
      const Scalar * a_data      = a.derived().data();
      const Scalar * b_data      = b.derived().data();
      const Scalar * c_data      = c.derived().data();
      Scalar *       cprime_data = m_cprime.data();

      cprime_data[0] = c_data[0] / b_data[0];

      // Loop ottimizzato con accesso sequenziale
      for ( integer i = 1; i < m_size - 1; ++i )
      {
        const Scalar denom = b_data[i] - a_data[i - 1] * cprime_data[i - 1];
        cprime_data[i]     = c_data[i] / denom;
      }
    }

    /*============================================================
      SOLVE - OTTIMIZZATA
    ============================================================*/

    // Versione template unica - nessun overload
    template <typename DerivedA, typename DerivedB, typename DerivedRHS, typename DerivedX> void solve(
      const Eigen::DenseBase<DerivedA> &   a,
      const Eigen::DenseBase<DerivedB> &   b,
      const Eigen::DenseBase<DerivedRHS> & rhs,
      Eigen::DenseBase<DerivedX> &         x )
    {
      assert(
        static_cast<integer>( a.size() ) == m_size - 1 && static_cast<integer>( b.size() ) == m_size &&
        static_cast<integer>( rhs.size() ) == m_size );

      if ( m_size == 0 ) return;

      x.derived().resize( m_size );

      if ( m_size == 1 )
      {
        x( 0 ) = rhs( 0 ) / b( 0 );
        return;
      }

      // Accesso diretto ai dati
      const Scalar * a_data      = a.derived().data();
      const Scalar * b_data      = b.derived().data();
      const Scalar * rhs_data    = rhs.derived().data();
      const Scalar * cprime_data = m_cprime.data();
      Scalar *       dprime_data = m_dprime.data();
      Scalar *       x_data      = x.derived().data();

      // Forward substitution
      dprime_data[0] = rhs_data[0] / b_data[0];
      for ( integer i = 1; i < m_size; ++i )
      {
        const Scalar denom = b_data[i] - a_data[i - 1] * cprime_data[i - 1];
        dprime_data[i]     = ( rhs_data[i] - a_data[i - 1] * dprime_data[i - 1] ) / denom;
      }

      // Backward substitution
      x_data[m_size - 1] = dprime_data[m_size - 1];
      for ( integer i = m_size - 2; i >= 0; --i ) { x_data[i] = dprime_data[i] - cprime_data[i] * x_data[i + 1]; }
    }

    // Versione che risolve in-place sovrascrivendo il RHS
    template <typename DerivedA, typename DerivedB, typename DerivedRHS> void solve_inplace(
      const Eigen::DenseBase<DerivedA> & a,
      const Eigen::DenseBase<DerivedB> & b,
      Eigen::DenseBase<DerivedRHS> &     rhs )
    {
      assert(
        static_cast<integer>( a.size() ) == m_size - 1 && static_cast<integer>( b.size() ) == m_size &&
        static_cast<integer>( rhs.size() ) == m_size );

      if ( m_size == 0 ) return;

      if ( m_size == 1 )
      {
        rhs( 0 ) /= b( 0 );
        return;
      }

      // Accesso diretto ai dati
      const Scalar * a_data      = a.derived().data();
      const Scalar * b_data      = b.derived().data();
      const Scalar * cprime_data = m_cprime.data();
      Scalar *       dprime_data = m_dprime.data();
      Scalar *       rhs_data    = rhs.derived().data();

      // Forward substitution - sovrascrive rhs con dprime
      rhs_data[0] /= b_data[0];
      dprime_data[0] = rhs_data[0];

      for ( integer i = 1; i < m_size; ++i )
      {
        const Scalar denom = b_data[i] - a_data[i - 1] * cprime_data[i - 1];
        rhs_data[i]        = ( rhs_data[i] - a_data[i - 1] * dprime_data[i - 1] ) / denom;
        dprime_data[i]     = rhs_data[i];
      }

      // Backward substitution - sovrascrive rhs con la soluzione
      for ( integer i = m_size - 2; i >= 0; --i ) { rhs_data[i] = dprime_data[i] - cprime_data[i] * rhs_data[i + 1]; }
    }

    // Versione che restituisce la soluzione
    template <typename DerivedA, typename DerivedB, typename DerivedRHS> VecS solve(
      const Eigen::DenseBase<DerivedA> &   a,
      const Eigen::DenseBase<DerivedB> &   b,
      const Eigen::DenseBase<DerivedRHS> & rhs )
    {
      VecS x;
      solve( a, b, rhs, x );
      return x;
    }

    // Batch solve ottimizzato con pre-allocazione
    template <typename DerivedA, typename DerivedB, typename DerivedRHS, typename DerivedX> void solve_batch(
      const Eigen::DenseBase<DerivedA> &   a,
      const Eigen::DenseBase<DerivedB> &   b,
      const Eigen::DenseBase<DerivedRHS> & RHS,
      Eigen::DenseBase<DerivedX> &         X )
    {
      const integer k = RHS.cols();
      X.derived().resize( m_size, k );

      // Pre-allocazione buffer riutilizzabile
      VecS x_col( m_size );

      // Loop interchange per migliore località cache
      for ( integer j = 0; j < k; ++j )
      {
        solve( a, b, RHS.col( j ), x_col );
        X.col( j ) = x_col;
      }
    }

    // Batch solve in-place
    template <typename DerivedA, typename DerivedB, typename DerivedX> void solve_batch_inplace(
      const Eigen::DenseBase<DerivedA> & a,
      const Eigen::DenseBase<DerivedB> & b,
      Eigen::DenseBase<DerivedX> &       X )
    {
      const integer k = X.cols();

      // Loop interchange per migliore località cache
      for ( integer j = 0; j < k; ++j )
      {
        auto col = X.col( j );
        solve_inplace( a, b, col );
      }
    }

    // Versione che restituisce la soluzione batch
    template <typename DerivedA, typename DerivedB, typename DerivedRHS> MatB solve_batch(
      const Eigen::DenseBase<DerivedA> &   a,
      const Eigen::DenseBase<DerivedB> &   b,
      const Eigen::DenseBase<DerivedRHS> & RHS )
    {
      MatB X;
      solve_batch( a, b, RHS, X );
      return X;
    }

    /*============================================================
      CYCLIC SOLVE - MANTIENE LOGICA ORIGINALE
    ============================================================*/

    // Versione template unica - nessun overload
    template <typename DerivedA, typename DerivedB, typename DerivedC, typename DerivedRHS, typename DerivedX>
    void solve_cyclic(
      const Eigen::DenseBase<DerivedA> &   a,
      const Eigen::DenseBase<DerivedB> &   b,
      const Eigen::DenseBase<DerivedC> &   c,
      Scalar                               alpha,
      Scalar                               beta,
      const Eigen::DenseBase<DerivedRHS> & rhs,
      Eigen::DenseBase<DerivedX> &         x )
    {
      x.derived().resize( m_size );

      if ( m_size < 3 )
      {
        if ( m_size == 1 )
        {
          x( 0 ) = rhs( 0 ) / ( b( 0 ) + alpha + beta );
          return;
        }
        else if ( m_size == 2 )
        {
          Scalar det = ( b( 0 ) + alpha ) * ( b( 1 ) + beta ) - c( 0 ) * a( 0 );
          x( 0 )     = ( ( b( 1 ) + beta ) * rhs( 0 ) - c( 0 ) * rhs( 1 ) ) / det;
          x( 1 )     = ( ( b( 0 ) + alpha ) * rhs( 1 ) - a( 0 ) * rhs( 0 ) ) / det;
          return;
        }
      }

      VecS bb = b;
      bb( 0 ) -= alpha;
      bb( m_size - 1 ) -= beta;

      factorize( a, bb, c );

      VecS uscalar( m_size );
      uscalar.setZero();
      uscalar( 0 )          = alpha;
      uscalar( m_size - 1 ) = beta;

      VecS Yscalar( m_size ), Zscalar( m_size );
      solve( a, bb, rhs, Yscalar );
      solve( a, bb, uscalar, Zscalar );

      Scalar fact = ( Yscalar( 0 ) + Yscalar( m_size - 1 ) ) / ( Scalar( 1 ) + Zscalar( 0 ) + Zscalar( m_size - 1 ) );
      x.derived() = Yscalar - fact * Zscalar;
    }

    // Versione che restituisce la soluzione
    template <typename DerivedA, typename DerivedB, typename DerivedC, typename DerivedRHS> VecS solve_cyclic(
      const Eigen::DenseBase<DerivedA> &   a,
      const Eigen::DenseBase<DerivedB> &   b,
      const Eigen::DenseBase<DerivedC> &   c,
      Scalar                               alpha,
      Scalar                               beta,
      const Eigen::DenseBase<DerivedRHS> & rhs )
    {
      VecS x;
      solve_cyclic( a, b, c, alpha, beta, rhs, x );
      return x;
    }
  };

  // ===========================================================================
  // SOLUTORE TRIDIAGONALE A BLOCCHI
  // ===========================================================================
  template <typename Scalar, int BlockSize = 1> class BlockTridiagonalSolver
  {
  public:
    static constexpr bool m_dynamic_block = BlockSize < 0;

    using Block = Eigen::Matrix<
      Scalar,
      ( m_dynamic_block ? Eigen::Dynamic : BlockSize ),
      ( m_dynamic_block ? Eigen::Dynamic : BlockSize )>;

    using VecB    = Eigen::Matrix<Scalar, ( m_dynamic_block ? Eigen::Dynamic : BlockSize ), 1>;
    using VecS    = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using MatB    = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using integer = Eigen::Index;

  private:
    integer m_size = 0;
    integer m_dim  = BlockSize > 0 ? BlockSize : 0;

    // Block workspace
    std::vector<Block> m_Cprime;
    std::vector<VecB>  m_Dprime;

    // Pre-allocazione decomposizioni LU per riutilizzo
    using LuDecomposition = Eigen::PartialPivLU<Block>;
    std::vector<LuDecomposition> m_lu_cache;
    bool                         m_lu_cache_valid = false;

    // Cyclic factorizations
    bool m_cyclic_1_factorized  = false;
    using Cyclic1LuType         = std::shared_ptr<Eigen::PartialPivLU<Block>>;
    Cyclic1LuType m_cyclic_1_lu = nullptr;

    bool m_cyclic_2_factorized  = false;
    using Cyclic2LuType         = std::shared_ptr<Eigen::PartialPivLU<MatB>>;
    Cyclic2LuType m_cyclic_2_lu = nullptr;
    MatB          m_cyclic_2_matrix;

    // Workspace temporaneo per evitare allocazioni ripetute
    mutable Block m_temp_block;
    mutable VecB  m_temp_vec;

    constexpr integer get_block_dim( integer m_ ) const noexcept
    {
      if constexpr ( m_dynamic_block )
        return m_;
      else
        return BlockSize;
    }

    // Versioni inline per operazioni critiche
    inline VecS block_vector_to_flat( const std::vector<VecB> & block_vec ) const
    {
      VecS     flat( m_size * m_dim );
      Scalar * data = flat.data();
      for ( integer i = 0; i < m_size; ++i )
      {
        const VecB & block = block_vec[i];
        std::memcpy( data + i * m_dim, block.data(), m_dim * sizeof( Scalar ) );
      }
      return flat;
    }

    inline std::vector<VecB> flat_to_block_vector( const VecS & flat ) const
    {
      std::vector<VecB> block_vec( m_size );
      const Scalar *    data = flat.data();
      for ( integer i = 0; i < m_size; ++i )
      {
        block_vec[i].resize( m_dim );
        std::memcpy( block_vec[i].data(), data + i * m_dim, m_dim * sizeof( Scalar ) );
      }
      return block_vec;
    }

  public:
    BlockTridiagonalSolver() noexcept {}

    BlockTridiagonalSolver( integer n_, integer m_ = BlockSize ) { resize( n_, m_ ); }

    void resize( integer n_, integer m_ = BlockSize )
    {
      m_size = n_;
      m_dim  = get_block_dim( m_ );

      m_cyclic_1_factorized = false;
      m_cyclic_2_factorized = false;
      m_lu_cache_valid      = false;

      m_cyclic_1_lu.reset();
      m_cyclic_2_lu.reset();

      // Caso a blocchi: allocate block arrays
      if ( m_size > 1 ) m_Cprime.resize( m_size - 1 );
      m_Dprime.resize( m_size );

      // Pre-allocazione cache LU
      m_lu_cache.clear();
      m_lu_cache.reserve( m_size );

      if constexpr ( m_dynamic_block )
      {
        for ( auto & c : m_Cprime ) c.resize( m_dim, m_dim );
        for ( auto & v : m_Dprime ) v.resize( m_dim );

        // Pre-allocazione workspace
        m_temp_block.resize( m_dim, m_dim );
        m_temp_vec.resize( m_dim );
      }
      else
      {
        // Per blocchi statici, le dimensioni sono già fissate
        m_temp_block = Block::Zero( m_dim, m_dim );
        m_temp_vec   = VecB::Zero( m_dim );
      }
    }

    integer size() const noexcept { return m_size; }
    integer block_dim() const noexcept { return m_dim; }

    /*============================================================
      FACTORIZATION - OTTIMIZZATA
    ============================================================*/

    // Riutilizzo decomposizioni LU quando possibile
    void factorize( const std::vector<Block> & A, const std::vector<Block> & B, const std::vector<Block> & C )
    {
      assert(
        static_cast<integer>( A.size() ) == m_size - 1 && static_cast<integer>( B.size() ) == m_size &&
        static_cast<integer>( C.size() ) == m_size - 1 );

      if ( m_size < 2 ) return;

      // Pre-allocazione LU cache se necessario
      if ( m_lu_cache.size() != static_cast<size_t>( m_size ) )
      {
        m_lu_cache.clear();
        m_lu_cache.reserve( m_size );
        for ( integer i = 0; i < m_size; ++i )
        {
          m_lu_cache.emplace_back( B[0] );  // Dummy initialization
        }
      }

      // Compute e memorizza LU per B[0]
      m_lu_cache[0].compute( B[0] );
      m_Cprime[0] = m_lu_cache[0].solve( C[0] );

      for ( integer i = 1; i < m_size - 1; ++i )
      {
        // Riutilizzo m_temp_block per evitare allocazioni
        m_temp_block = B[i];
        m_temp_block.noalias() -= A[i - 1] * m_Cprime[i - 1];

        m_lu_cache[i].compute( m_temp_block );
        m_Cprime[i] = m_lu_cache[i].solve( C[i] );
      }

      m_lu_cache_valid = true;
    }

    /*============================================================
      SOLVE - OTTIMIZZATA
    ============================================================*/

    // Utilizzo LU cache pre-calcolata
    void solve(
      const std::vector<Block> & A,
      const std::vector<Block> & B,
      const std::vector<VecB> &  RHS,
      std::vector<VecB> &        X )
    {
      assert(
        static_cast<integer>( A.size() ) == m_size - 1 && static_cast<integer>( B.size() ) == m_size &&
        static_cast<integer>( RHS.size() ) == m_size );

      if ( m_size == 0 ) return;

      X.resize( m_size );

      if ( m_size == 1 )
      {
        if ( m_lu_cache_valid && !m_lu_cache.empty() ) { X[0] = m_lu_cache[0].solve( RHS[0] ); }
        else
        {
          X[0] = B[0].partialPivLu().solve( RHS[0] );
        }
        return;
      }

      // Forward substitution con riutilizzo LU cache
      if ( m_lu_cache_valid ) { m_Dprime[0] = m_lu_cache[0].solve( RHS[0] ); }
      else
      {
        m_Dprime[0] = B[0].partialPivLu().solve( RHS[0] );
      }

      for ( integer i = 1; i < m_size; ++i )
      {
        m_temp_vec = RHS[i];
        m_temp_vec.noalias() -= A[i - 1] * m_Dprime[i - 1];

        if ( m_lu_cache_valid && i < m_size - 1 ) { m_Dprime[i] = m_lu_cache[i].solve( m_temp_vec ); }
        else
        {
          m_temp_block = B[i];
          m_temp_block.noalias() -= A[i - 1] * m_Cprime[i - 1];
          m_Dprime[i] = m_temp_block.partialPivLu().solve( m_temp_vec );
        }
      }

      // Backward substitution
      X[m_size - 1] = m_Dprime[m_size - 1];
      for ( integer i = m_size - 2; i >= 0; --i )
      {
        X[i] = m_Dprime[i];
        X[i].noalias() -= m_Cprime[i] * X[i + 1];
      }
    }

    // Versione che risolve in-place sovrascrivendo RHS
    void solve_inplace(
      const std::vector<Block> & A,
      const std::vector<Block> & B,
      std::vector<VecB> &        RHS_X )  // Input: RHS, Output: soluzione
    {
      assert(
        static_cast<integer>( A.size() ) == m_size - 1 && static_cast<integer>( B.size() ) == m_size &&
        static_cast<integer>( RHS_X.size() ) == m_size );

      if ( m_size == 0 ) return;

      if ( m_size == 1 )
      {
        if ( m_lu_cache_valid && !m_lu_cache.empty() ) { RHS_X[0] = m_lu_cache[0].solve( RHS_X[0] ); }
        else
        {
          RHS_X[0] = B[0].partialPivLu().solve( RHS_X[0] );
        }
        return;
      }

      // Forward substitution con riutilizzo LU cache
      // Usiamo RHS_X come workspace per Dprime
      if ( m_lu_cache_valid ) { m_Dprime[0] = m_lu_cache[0].solve( RHS_X[0] ); }
      else
      {
        m_Dprime[0] = B[0].partialPivLu().solve( RHS_X[0] );
      }
      RHS_X[0] = m_Dprime[0];  // Sovrascrivi con il valore di Dprime

      for ( integer i = 1; i < m_size; ++i )
      {
        m_temp_vec = RHS_X[i];
        m_temp_vec.noalias() -= A[i - 1] * m_Dprime[i - 1];

        if ( m_lu_cache_valid && i < m_size - 1 ) { m_Dprime[i] = m_lu_cache[i].solve( m_temp_vec ); }
        else
        {
          m_temp_block = B[i];
          m_temp_block.noalias() -= A[i - 1] * m_Cprime[i - 1];
          m_Dprime[i] = m_temp_block.partialPivLu().solve( m_temp_vec );
        }
        RHS_X[i] = m_Dprime[i];  // Sovrascrivi con il valore di Dprime
      }

      // Backward substitution - sovrascrive RHS_X con la soluzione
      // RHS_X[m_size-1] già contiene Dprime[m_size-1]
      for ( integer i = m_size - 2; i >= 0; --i )
      {
        RHS_X[i] = m_Dprime[i];
        RHS_X[i].noalias() -= m_Cprime[i] * RHS_X[i + 1];
      }
    }

    // Versione che restituisce la soluzione
    std::vector<VecB> solve( const std::vector<Block> & A, const std::vector<Block> & B, const std::vector<VecB> & RHS )
    {
      std::vector<VecB> X;
      solve( A, B, RHS, X );
      return X;
    }

    /*============================================================
      CYCLIC SOLVE
    ============================================================*/

    void factorize_cyclic_special(
      const std::vector<Block> & A,
      const std::vector<Block> & B,
      const std::vector<Block> & C,
      const Block &              Alpha,
      const Block &              Beta )
    {
      m_cyclic_1_factorized = false;
      m_cyclic_2_factorized = false;
      m_cyclic_1_lu.reset();
      m_cyclic_2_lu.reset();

      if ( m_size == 1 )
      {
        Block Btotal          = B[0] + Alpha + Beta;
        m_cyclic_1_lu         = std::make_shared<Eigen::PartialPivLU<Block>>( Btotal );
        m_cyclic_1_factorized = true;
        return;
      }

      if ( m_size == 2 )
      {
        integer m = B[0].rows();
        m_cyclic_2_matrix.resize( 2 * m, 2 * m );
        m_cyclic_2_matrix.setZero();

        m_cyclic_2_matrix.block( 0, 0, m, m ) = B[0];
        if ( C.empty() ) { m_cyclic_2_matrix.block( 0, m, m, m ) = Alpha; }
        else
        {
          m_cyclic_2_matrix.block( 0, m, m, m ) = C[0];
        }
        if ( A.empty() ) { m_cyclic_2_matrix.block( m, 0, m, m ) = Beta; }
        else
        {
          m_cyclic_2_matrix.block( m, 0, m, m ) = A[0];
        }
        m_cyclic_2_matrix.block( m, m, m, m ) = B[1];

        m_cyclic_2_lu         = std::make_shared<Eigen::PartialPivLU<MatB>>( m_cyclic_2_matrix );
        m_cyclic_2_factorized = true;
        return;
      }
    }

    void solve_cyclic_using_factorization( const std::vector<VecB> & RHS, std::vector<VecB> & X )
    {
      assert( static_cast<integer>( RHS.size() ) == m_size );

      X.resize( m_size );

      if ( m_size == 1 )
      {
        if ( !m_cyclic_1_factorized || !m_cyclic_1_lu )
          throw std::runtime_error( "Cyclic factorization for n=1 not available" );
        X[0] = m_cyclic_1_lu->solve( RHS[0] );
        return;
      }

      if ( m_size == 2 )
      {
        if ( !m_cyclic_2_factorized || !m_cyclic_2_lu )
          throw std::runtime_error( "Cyclic factorization for n=2 not available" );

        integer m = RHS[0].size();
        VecS    rhs_full( 2 * m );
        rhs_full.segment( 0, m ) = RHS[0];
        rhs_full.segment( m, m ) = RHS[1];

        VecS sol_full = m_cyclic_2_lu->solve( rhs_full );

        X[0] = sol_full.segment( 0, m );
        X[1] = sol_full.segment( m, m );
        return;
      }

      throw std::runtime_error( "solve_cyclic_using_factorization() only supports n=1 and n=2" );
    }

    void solve_cyclic(
      const std::vector<Block> & A,
      const std::vector<Block> & B,
      const std::vector<Block> & C,
      const Block &              Alpha,
      const Block &              Beta,
      const std::vector<VecB> &  RHS,
      std::vector<VecB> &        X )
    {
      assert( static_cast<integer>( RHS.size() ) == m_size );

      if ( m_size == 1 || m_size == 2 )
      {
        if ( ( m_size == 1 && m_cyclic_1_factorized ) || ( m_size == 2 && m_cyclic_2_factorized ) )
        {
          solve_cyclic_using_factorization( RHS, X );
          return;
        }

        if ( m_size == 1 )
        {
          X.resize( 1 );
          Block Btotal = B[0] + Alpha + Beta;
          X[0]         = Btotal.partialPivLu().solve( RHS[0] );
          return;
        }

        if ( m_size == 2 )
        {
          X.resize( 2 );
          integer m = B[0].rows();
          MatB    M_sys( 2 * m, 2 * m );
          M_sys.setZero();

          M_sys.block( 0, 0, m, m ) = B[0];
          M_sys.block( m, m, m, m ) = B[1];
          if ( C.empty() ) { M_sys.block( 0, m, m, m ) = Alpha; }
          else
          {
            M_sys.block( 0, m, m, m ) = C[0];
          }
          if ( A.empty() ) { M_sys.block( m, 0, m, m ) = Beta; }
          else
          {
            M_sys.block( m, 0, m, m ) = A[0];
          }

          VecS R_sys( 2 * m );
          R_sys.segment( 0, m ) = RHS[0];
          R_sys.segment( m, m ) = RHS[1];

          VecS Sol = M_sys.partialPivLu().solve( R_sys );
          X[0]     = Sol.segment( 0, m );
          X[1]     = Sol.segment( m, m );
          return;
        }
      }

      // General case
      std::vector<Block> BB = B;
      factorize( A, BB, C );

      std::vector<VecB> Y;
      solve( A, BB, RHS, Y );

      const integer m     = m_dim;
      const integer two_m = 2 * m;

      MatB U_mat( m_size * m, two_m ), V_mat( m_size * m, two_m );
      U_mat.setZero();
      V_mat.setZero();

      U_mat.block( 0, 0, m, m )                  = Block::Identity( m, m );
      U_mat.block( ( m_size - 1 ) * m, m, m, m ) = Block::Identity( m, m );
      V_mat.block( ( m_size - 1 ) * m, 0, m, m ) = Alpha.transpose();
      V_mat.block( 0, m, m, m )                  = Beta.transpose();

      MatB Z_mat( m_size * m, two_m );

      for ( integer col = 0; col < two_m; ++col )
      {
        std::vector<VecB> u_col( m_size ), z_col( m_size );
        for ( integer i = 0; i < m_size; ++i ) { u_col[i] = U_mat.block( i * m, col, m, 1 ); }

        solve( A, BB, u_col, z_col );

        for ( integer i = 0; i < m_size; ++i ) { Z_mat.block( i * m, col, m, 1 ) = z_col[i]; }
      }

      VecS Y_vec    = block_vector_to_flat( Y );
      VecS VtY      = V_mat.transpose() * Y_vec;
      MatB VtZ      = V_mat.transpose() * Z_mat;
      MatB IplusVtZ = MatB::Identity( two_m, two_m ) + VtZ;
      VecS K        = IplusVtZ.partialPivLu().solve( VtY );
      VecS X_vec    = Y_vec - Z_mat * K;

      X = flat_to_block_vector( X_vec );
    }

    // Versione che restituisce la soluzione
    std::vector<VecB> solve_cyclic(
      const std::vector<Block> & A,
      const std::vector<Block> & B,
      const std::vector<Block> & C,
      const Block &              Alpha,
      const Block &              Beta,
      const std::vector<VecB> &  RHS )
    {
      std::vector<VecB> X;
      solve_cyclic( A, B, C, Alpha, Beta, RHS, X );
      return X;
    }
  };

}  // namespace Utils

#endif

// EOF: Utils_TridiagonalSolver.hh
