/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2025                                                      |
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
 |      Università degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

//
// file: Utils_ssolver.hh
//

#pragma once

#ifndef UTILS_SSOLVER_dot_HH
#define UTILS_SSOLVER_dot_HH

#include "Utils.hh"
#include "Utils_eigen.hh"
#include "Utils_fmt.hh"

namespace Utils
{

  using std::abs;
  using std::max;
  using std::min;
  using std::pow;
  using std::sqrt;

  /**
   * @class DenseSymmetricSolver
   * @brief Solver for dense symmetric linear systems of the form
   * \f$(\bm{A} + \lambda \bm{I})\bm{x} = \bm{b}\f$.
   */
  template <typename Scalar = double> class DenseSymmetricSolver
  {
  private:
    using integer      = Eigen::Index;
    using Matrix       = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using SparseMatrix = Eigen::SparseMatrix<Scalar>;
    using Vector       = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using Solver       = Eigen::LDLT<Matrix>;

    integer m_dim;
    Solver  m_solver;

  public:
    explicit DenseSymmetricSolver( Matrix const & A, Scalar lambda )
    {
      UTILS_ASSERT( A.rows() == A.cols(), "DenseSymmetricSolver( A, λ={} ) Matrix A must be square.", lambda );
      UTILS_ASSERT( Utils::isSymmetric( A ), "DenseSymmetricSolver( A, λ={} ) Matrix A must be symmetric.", lambda );
      m_dim = A.rows();
      new ( &m_solver ) Solver( A + lambda * Matrix::Identity( A.rows(), A.cols() ) );
    }

    explicit DenseSymmetricSolver( SparseMatrix const & A, Scalar lambda )
    {
      UTILS_ASSERT( A.rows() == A.cols(), "DenseSymmetricSolver( A, λ={} ) Matrix A must be square.", lambda );
      UTILS_ASSERT( Utils::isSymmetric( A ), "DenseSymmetricSolver( A, λ={} ) Matrix A must be symmetric.", lambda );
      m_dim     = A.rows();
      Matrix AA = A;
      AA += lambda * Matrix::Identity( A.rows(), A.cols() );
      new ( &m_solver ) Solver( AA );
    }

    Vector solve( Vector const & b ) const
    {
      UTILS_ASSERT(
        b.size() == m_dim,
        "DenseSymmetricSolver.solve( b ) b has incorrect dimension #b = {} expected {}",
        b.size(),
        m_dim );
      return m_solver.solve( b );
    }

    void solve_in_place( Vector & b ) const
    {
      UTILS_ASSERT(
        b.size() == m_dim,
        "DenseSymmetricSolver.solve_in_place( b ) b has incorrect dimension #b = {} expected {}",
        b.size(),
        m_dim );
      b = m_solver.solve( b );
    }
  };

  /**
   * @class SparseSymmetricSolver
   * @brief Solver for sparse symmetric linear systems of the form
   * \f$(\bm{A} + \lambda \bm{I})\bm{x} = \bm{b}\f$.
   */
  template <typename Scalar = double> class SparseSymmetricSolver
  {
  private:
    using integer      = Eigen::Index;
    using Matrix       = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using SparseMatrix = Eigen::SparseMatrix<Scalar>;
    using Vector       = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using Solver       = Eigen::SimplicialLDLT<SparseMatrix>;

    integer m_dim;
    Solver  m_solver;

  public:
    explicit SparseSymmetricSolver( SparseMatrix const & A, Scalar lambda )
    {
      UTILS_ASSERT( A.rows() == A.cols(), "SparseSymmetricSolver( A, λ={} ) Matrix A must be square.", lambda );
      UTILS_ASSERT( Utils::isSymmetric( A ), "SparseSymmetricSolver( A, λ={} ) Matrix A must be symmetric.", lambda );
      m_dim = A.rows();

      SparseMatrix M = A;
      for ( integer i = 0; i < m_dim; ++i ) M.coeffRef( i, i ) += lambda;
      M.makeCompressed();

      m_solver.compute( M );

      UTILS_ASSERT(
        m_solver.info() == Eigen::Success,
        "SparseSymmetricSolver( A, λ={} ) Factorization failed. Matrix may be singular or indefinite.",
        lambda );
    }

    explicit SparseSymmetricSolver( Matrix const & A, Scalar lambda )
    {
      UTILS_ASSERT( A.rows() == A.cols(), "SparseSymmetricSolver( A, λ={} ) Matrix A must be square.", lambda );
      UTILS_ASSERT( Utils::isSymmetric( A ), "SparseSymmetricSolver( A, λ={} ) Matrix A must be symmetric.", lambda );
      m_dim = A.rows();

      SparseMatrix M = A;
      for ( integer i = 0; i < m_dim; ++i ) M.coeffRef( i, i ) += lambda;
      M.makeCompressed();

      m_solver.compute( M );

      UTILS_ASSERT(
        m_solver.info() == Eigen::Success,
        "SparseSymmetricSolver( A, λ={} ) Factorization failed. Matrix may be singular or indefinite.",
        lambda );
    }

    Vector solve( Vector const & b ) const
    {
      UTILS_ASSERT(
        b.size() == m_dim,
        "SparseSymmetricSolver.solve( b ) b has incorrect dimension #b = {} expected {}",
        b.size(),
        m_dim );
      return m_solver.solve( b );
    }

    void solve_in_place( Vector & b ) const
    {
      UTILS_ASSERT(
        b.size() == m_dim,
        "SparseSymmetricSolver.solve_in_place( b ) b has incorrect dimension #b = {} expected {}",
        b.size(),
        m_dim );
      b = m_solver.solve( b );
    }
  };

  /**
   * @class SparsePatternAnalysis
   * @brief Analyzes sparse matrix pattern to detect special structures
   */
  template <typename Scalar = double> class SparsePatternAnalysis
  {
  private:
    using integer      = Eigen::Index;
    using SparseMatrix = Eigen::SparseMatrix<Scalar>;

  public:
    struct DiagonalBlock
    {
      integer start_row;
      integer size;
      bool    is_dense{ true };

      DiagonalBlock( integer sr, integer sz, bool dense = true ) : start_row( sr ), size( sz ), is_dense( dense ) {}
    };

    struct Partitioning
    {
      bool                       is_purely_diagonal{ false };
      std::vector<DiagonalBlock> diagonal_blocks;

      bool has_partitioning() const { return is_purely_diagonal || !diagonal_blocks.empty(); }
    };

    static bool is_purely_diagonal( const SparseMatrix & A )
    {
      for ( int k = 0; k < A.outerSize(); ++k )
      {
        for ( typename SparseMatrix::InnerIterator it( A, k ); it; ++it )
        {
          if ( it.row() != it.col() ) return false;
        }
      }
      return true;
    }

    static Partitioning analyze_pattern( const SparseMatrix & A, integer dense_threshold = 100 )
    {
      Partitioning partitioning;
      integer      n = A.rows();

      if ( is_purely_diagonal( A ) )
      {
        partitioning.is_purely_diagonal = true;
        return partitioning;
      }

      std::vector<std::vector<integer>> adj( n );
      for ( int k = 0; k < A.outerSize(); ++k )
      {
        for ( typename SparseMatrix::InnerIterator it( A, k ); it; ++it )
        {
          integer i = it.row();
          integer j = it.col();
          if ( i != j && it.value() != 0 )
          {
            adj[i].push_back( j );
            adj[j].push_back( i );
          }
        }
      }

      std::vector<bool> visited( n, false );

      for ( integer i = 0; i < n; ++i )
      {
        if ( !visited[i] )
        {
          std::vector<integer> block_indices;
          std::queue<integer>  q;
          q.push( i );
          visited[i] = true;

          while ( !q.empty() )
          {
            integer curr = q.front();
            q.pop();
            block_indices.push_back( curr );

            for ( integer neighbor : adj[curr] )
            {
              if ( !visited[neighbor] )
              {
                visited[neighbor] = true;
                q.push( neighbor );
              }
            }
          }

          std::sort( block_indices.begin(), block_indices.end() );

          bool is_contiguous = true;
          for ( size_t k = 1; k < block_indices.size(); ++k )
          {
            if ( block_indices[k] != block_indices[k - 1] + 1 )
            {
              is_contiguous = false;
              break;
            }
          }

          if ( is_contiguous && !block_indices.empty() )
          {
            DiagonalBlock block(
              block_indices.front(),
              static_cast<integer>( block_indices.size() ),
              static_cast<integer>( block_indices.size() ) <= dense_threshold );
            partitioning.diagonal_blocks.push_back( block );
          }
        }
      }

      if (
        partitioning.diagonal_blocks.size() == 1 && partitioning.diagonal_blocks[0].start_row == 0 &&
        partitioning.diagonal_blocks[0].size == n )
      {
        partitioning.diagonal_blocks.clear();
      }

      return partitioning;
    }
  };

  /**
   * @class SymmetricSolver
   * @brief Adaptive solver for symmetric linear systems
   */
  template <typename Scalar = double, Eigen::Index DenseThreshold = 100> class SymmetricSolver
  {
  private:
    using integer       = Eigen::Index;
    using Matrix        = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using SparseMatrix  = Eigen::SparseMatrix<Scalar>;
    using Vector        = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using Partitioning  = typename SparsePatternAnalysis<Scalar>::Partitioning;
    using DiagonalBlock = typename SparsePatternAnalysis<Scalar>::DiagonalBlock;

    using DenseSolver  = DenseSymmetricSolver<Scalar>;
    using SparseSolver = SparseSymmetricSolver<Scalar>;

    struct BlockSolver
    {
      DiagonalBlock                 block;
      std::unique_ptr<DenseSolver>  dense_solver;
      std::unique_ptr<SparseSolver> sparse_solver;

      BlockSolver( const DiagonalBlock & blk ) : block( blk ) {}

      BlockSolver( const BlockSolver & )                 = delete;
      BlockSolver & operator=( const BlockSolver & )     = delete;
      BlockSolver( BlockSolver && ) noexcept             = default;
      BlockSolver & operator=( BlockSolver && ) noexcept = default;

      Vector solve( const Vector & b_block ) const
      {
        if ( block.is_dense && dense_solver ) { return dense_solver->solve( b_block ); }
        else if ( !block.is_dense && sparse_solver ) { return sparse_solver->solve( b_block ); }
        UTILS_ASSERT( false, "BlockSolver: solver not initialized" );
        return Vector();
      }

      void solve_in_place( Vector & b_block ) const
      {
        UTILS_ASSERT(
          b_block.size() == block.size,
          "BlockSolver.solve_in_place: wrong dimension {} expected {}",
          b_block.size(),
          block.size );

        if ( block.is_dense && dense_solver ) { dense_solver->solve_in_place( b_block ); }
        else if ( !block.is_dense && sparse_solver ) { sparse_solver->solve_in_place( b_block ); }
        else
        {
          UTILS_ASSERT( false, "BlockSolver: solver not initialized" );
        }
      }
    };

    integer m_dim;
    Scalar  m_lambda;

    enum class Strategy
    {
      Dense,
      Sparse,
      Diagonal,
      BlockDiagonal
    };

    Strategy m_strategy;

    std::optional<DenseSolver>  m_dense_solver;
    std::optional<SparseSolver> m_sparse_solver;

    Vector m_diagonal_values;

    std::vector<BlockSolver> m_block_solvers;

    void initialize_from_dense( const Matrix & A, Scalar lambda )
    {
      m_dim    = A.rows();
      m_lambda = lambda;

      UTILS_ASSERT( A.rows() == A.cols(), "SymmetricSolver( A, λ={} ) Matrix A must be square.", lambda );
      UTILS_ASSERT( Utils::isSymmetric( A ), "SymmetricSolver( A, λ={} ) Matrix A must be symmetric.", lambda );

      if ( m_dim <= DenseThreshold )
      {
        m_strategy = Strategy::Dense;
        m_dense_solver.emplace( A, lambda );
      }
      else
      {
        m_strategy = Strategy::Sparse;
        m_sparse_solver.emplace( A, lambda );
      }
    }

    void initialize_from_sparse(
      const SparseMatrix & A,
      Scalar               lambda,
      const Partitioning & partitioning = Partitioning() )
    {
      m_dim    = A.rows();
      m_lambda = lambda;

      UTILS_ASSERT( A.rows() == A.cols(), "SymmetricSolver( A, λ={} ) Matrix A must be square.", lambda );
      UTILS_ASSERT( Utils::isSymmetric( A ), "SymmetricSolver( A, λ={} ) Matrix A must be symmetric.", lambda );

      Partitioning part = partitioning.has_partitioning()
                            ? partitioning
                            : SparsePatternAnalysis<Scalar>::analyze_pattern( A, DenseThreshold );

      if ( part.is_purely_diagonal )
      {
        m_strategy = Strategy::Diagonal;
        m_diagonal_values.resize( m_dim );

        for ( integer i = 0; i < m_dim; ++i ) { m_diagonal_values( i ) = A.coeff( i, i ) + lambda; }
        return;
      }

      if ( !part.diagonal_blocks.empty() )
      {
        m_strategy = Strategy::BlockDiagonal;

        // Verifica che i blocchi coprano l'intera matrice
        integer total_block_size = 0;
        for ( const auto & block : part.diagonal_blocks ) { total_block_size += block.size; }

        // Se i blocchi non coprono l'intera matrice, usa la strategia sparse
        if ( total_block_size != m_dim )
        {
          m_strategy = Strategy::Sparse;
          m_sparse_solver.emplace( A, lambda );
          return;
        }

        for ( const auto & block : part.diagonal_blocks )
        {
          BlockSolver block_solver( block );

          SparseMatrix subA = A.middleRows( block.start_row, block.size ).middleCols( block.start_row, block.size );

          if ( block.is_dense )
          {
            Matrix denseA             = subA;
            block_solver.dense_solver = std::make_unique<DenseSolver>( denseA, lambda );
          }
          else
          {
            block_solver.sparse_solver = std::make_unique<SparseSolver>( subA, lambda );
          }

          m_block_solvers.push_back( std::move( block_solver ) );
        }
        return;
      }

      // General sparse matrix - always use sparse solver
      m_strategy = Strategy::Sparse;
      m_sparse_solver.emplace( A, lambda );
    }

  public:
    explicit SymmetricSolver( const Matrix & A, Scalar lambda ) { initialize_from_dense( A, lambda ); }

    explicit SymmetricSolver(
      const SparseMatrix & A,
      Scalar               lambda,
      const Partitioning & partitioning = Partitioning() )
    {
      initialize_from_sparse( A, lambda, partitioning );
    }

    Vector solve( const Vector & b ) const
    {
      UTILS_ASSERT(
        b.size() == m_dim,
        "SymmetricSolver.solve(b): b has incorrect dimension #b = {} expected {}",
        b.size(),
        m_dim );

      switch ( m_strategy )
      {
        case Strategy::Dense: return m_dense_solver->solve( b );

        case Strategy::Sparse: return m_sparse_solver->solve( b );

        case Strategy::Diagonal:
        {
          Vector x( m_dim );
          for ( integer i = 0; i < m_dim; ++i ) { x( i ) = b( i ) / m_diagonal_values( i ); }
          return x;
        }

        case Strategy::BlockDiagonal:
        {
          Vector x = Vector::Zero( m_dim );

          for ( const auto & block_solver : m_block_solvers )
          {
            Vector b_block = b.segment( block_solver.block.start_row, block_solver.block.size );
            Vector x_block = block_solver.solve( b_block );
            x.segment( block_solver.block.start_row, block_solver.block.size ) = x_block;
          }
          return x;
        }
      }

      UTILS_ASSERT( false, "SymmetricSolver: invalid strategy" );
      return Vector();
    }

    void solve_in_place( Vector & b ) const
    {
      UTILS_ASSERT(
        b.size() == m_dim,
        "SymmetricSolver.solve_in_place(b): b has incorrect dimension #b = {} expected {}",
        b.size(),
        m_dim );

      switch ( m_strategy )
      {
        case Strategy::Dense: m_dense_solver->solve_in_place( b ); break;

        case Strategy::Sparse: m_sparse_solver->solve_in_place( b ); break;

        case Strategy::Diagonal:
          for ( integer i = 0; i < m_dim; ++i ) { b( i ) /= m_diagonal_values( i ); }
          break;

        case Strategy::BlockDiagonal:
          // Risolve ogni blocco direttamente sul segmento corrispondente
          for ( const auto & block_solver : m_block_solvers )
          {
            auto b_segment = b.segment( block_solver.block.start_row, block_solver.block.size );
            block_solver.solve_in_place( b_segment );
          }
          break;

        default: UTILS_ASSERT( false, "SymmetricSolver: invalid strategy" );
      }
    }

    std::string get_strategy_name() const
    {
      switch ( m_strategy )
      {
        case Strategy::Dense: return "Dense";
        case Strategy::Sparse: return "Sparse";
        case Strategy::Diagonal: return "Diagonal";
        case Strategy::BlockDiagonal: return "BlockDiagonal (" + std::to_string( m_block_solvers.size() ) + " blocks)";
        default: return "Unknown";
      }
    }

    bool is_block_diagonal() const { return m_strategy == Strategy::BlockDiagonal; }

    integer get_num_blocks() const { return static_cast<integer>( m_block_solvers.size() ); }

    // Metodo per parallelizzare la risoluzione dei blocchi (opzionale)
    template <bool UseParallel = false> Vector solve_parallel( const Vector & b ) const
    {
      static_assert( !UseParallel || DenseThreshold > 0, "Parallel execution requires DenseThreshold > 0" );

      if constexpr ( !UseParallel || m_strategy != Strategy::BlockDiagonal ) { return solve( b ); }
      else
      {
        UTILS_ASSERT(
          b.size() == m_dim,
          "SymmetricSolver.solve_parallel(b): b has incorrect dimension #b = {} expected {}",
          b.size(),
          m_dim );

        Vector x = Vector::Zero( m_dim );

        // Parallelizzazione semplice (da implementare con OpenMP o altro se necessario)
        for ( size_t i = 0; i < m_block_solvers.size(); ++i )
        {
          const auto & block_solver = m_block_solvers[i];
          Vector       b_block      = b.segment( block_solver.block.start_row, block_solver.block.size );
          Vector       x_block      = block_solver.solve( b_block );
          x.segment( block_solver.block.start_row, block_solver.block.size ) = x_block;
        }

        return x;
      }
    }
  };
}  // namespace Utils

#endif
