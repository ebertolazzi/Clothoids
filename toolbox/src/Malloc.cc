/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2020                                                      |
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
//
// file: Malloc.cc
//

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#if defined( __llvm__ ) || defined( __clang__ )
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#pragma clang diagnostic ignored "-Wduplicate-enum"
#endif

#include <iostream>

#include "Utils.hh"
#include "Utils_fmt.hh"
#include "Utils_trace.hh"

namespace Utils
{

  using std::cerr;
  using std::exception;
  using std::exit;
  using std::lock_guard;
  using std::mutex;

  std::mutex MallocMutex;

  int64_t CountAlloc{ 0 };
  int64_t CountFreed{ 0 };
  int64_t AllocatedBytes{ 0 };
  int64_t MaximumAllocatedBytes{ 0 };
  bool    MallocDebug{ false };

  string
  out_bytes( size_t nb )
  {
    size_t const Kb{ nb >> 10 };
    size_t const Mb{ Kb >> 10 };
    if ( size_t const Gb{ Mb >> 10 }; Gb > 0 )
    {
      size_t const mb{ ( 100 * ( Mb & 0x3FF ) ) / 1024 };
      return fmt::format( "{}Gb(+{}Mb)", Gb, mb );
    }
    if ( Mb > 0 )
    {
      size_t const kb{ ( 100 * ( Kb & 0x3FF ) ) / 1024 };
      return fmt::format( "{}Mb(+{}Kb)", Mb, kb );
    }
    if ( Kb > 0 )
    {
      size_t const b{ ( 100 * ( nb & 0x3FF ) ) / 1024 };
      return fmt::format( "{}Kb(+{}b)", Kb, b );
    }
    return fmt::format( "{}bytes", nb );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  Malloc<T>::allocate_internal( size_t const n )
  {
    try
    {
      size_t nb;
      {
        lock_guard lock( Utils::MallocMutex );
        nb = m_num_total_reserved * sizeof( T );
        ++CountFreed;
        AllocatedBytes -= nb;
      }

      delete[] m_p_memory;
      m_num_total_values   = n;
      m_num_total_reserved = n + ( n >> 3 );  // 12% more values
      m_p_memory           = new T[m_num_total_reserved];

      {
        lock_guard lock( Utils::MallocMutex );
        ++CountAlloc;
        nb = m_num_total_reserved * sizeof( T );
        AllocatedBytes += nb;
        if ( MaximumAllocatedBytes < AllocatedBytes ) MaximumAllocatedBytes = AllocatedBytes;
      }

      if ( MallocDebug ) fmt::print( "Allocating {} for {}\n", out_bytes( nb ), m_name );
    }
    catch ( exception const & exc )
    {
      string const reason = fmt::format( "Memory allocation failed: {}\nTry to allocate {} bytes for {}\n", exc.what(),
                                         n, m_name );
      print_trace( __LINE__, __FILE__, reason, cerr );
      exit( 0 );
    }
    catch ( ... )
    {
      string const reason = fmt::format( "Memory allocation failed for {}: memory exausted\n", m_name );
      print_trace( __LINE__, __FILE__, reason, cerr );
      exit( 0 );
    }
    m_num_total_values = n;
    m_num_allocated    = 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  Malloc<T>::allocate( size_t n )
  {
    UTILS_ASSERT( m_num_allocated == 0,
                  "Malloc[{}]::allocate( {} ), try to allocate already "
                  "allocated memory!\n",
                  m_name, n );
    if ( n > m_num_total_reserved ) allocate_internal( n );
    m_num_total_values = n;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  Malloc<T>::reallocate( size_t const n )
  {
    if ( n > m_num_total_reserved ) allocate_internal( n );
    m_num_total_values = n;
    m_num_allocated    = 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  T *
  Malloc<T>::malloc( size_t n )
  {
    UTILS_ASSERT( m_num_allocated == 0, "Malloc[{}]::malloc( {} ), try to allocate already allocated memory!\n", m_name,
                  n );
    if ( n > m_num_total_reserved ) allocate_internal( n );
    m_num_total_values = n;
    m_num_allocated    = n;
    return m_p_memory;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  T *
  Malloc<T>::realloc( size_t const n )
  {
    if ( n > m_num_total_reserved ) allocate_internal( n );
    m_num_total_values = n;
    m_num_allocated    = n;
    return m_p_memory;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  Malloc<T>::hard_free()
  {
    if ( m_p_memory != nullptr )
    {
      size_t nb;
      {
        lock_guard lock( Utils::MallocMutex );
        nb = m_num_total_reserved * sizeof( T );
        ++CountFreed;
        AllocatedBytes -= nb;
      }

      if ( MallocDebug ) fmt::print( "Freeing {} for {}\n", out_bytes( nb ), m_name );

      delete[] m_p_memory;
      m_p_memory           = nullptr;
      m_num_total_values   = 0;
      m_num_total_reserved = 0;
      m_num_allocated      = 0;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  Malloc<T>::memory_exausted( size_t sz )
  {
    string const reason = fmt::format( "Malloc<{}>::operator () ({}) -- Memory EXAUSTED\n", m_name, sz );
    print_trace( __LINE__, __FILE__, reason, cerr );
    exit( 0 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  Malloc<T>::pop_exausted( size_t sz )
  {
    string const reason = fmt::format( "Malloc<{}>::pop({}) -- Not enough element on Stack\n", m_name, sz );
    print_trace( __LINE__, __FILE__, reason, cerr );
    exit( 0 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  Malloc<T>::must_be_empty( string_view const where ) const
  {
    if ( m_num_allocated < m_num_total_values )
    {
      string const tmp = fmt::format( "in {} {}: not fully used!\nUnused: {} values\n", m_name, where,
                                      m_num_total_values - m_num_allocated );
      print_trace( __LINE__, __FILE__, tmp, cerr );
    }
    if ( m_num_allocated > m_num_total_values )
    {
      string const tmp = fmt::format( "in {} {}: too much used!\nMore used: {} values\n", m_name, where,
                                      m_num_allocated - m_num_total_values );
      print_trace( __LINE__, __FILE__, tmp, cerr );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  std::string
  Malloc<T>::info( string_view where ) const
  {
    std::size_t diff = m_num_allocated > m_num_total_values ? m_num_allocated - m_num_total_values
                                                            : m_num_total_values - m_num_allocated;
    return fmt::format(
        "in {} {}\n"
        "Allocated:  {}\n"
        "Reserved:   {}\n"
        "Difference: {} [|A-R|]\n",
        m_name, where, m_num_allocated, m_num_total_values, diff );
  }

  template class Malloc<char>;
  template class Malloc<uint16_t>;
  template class Malloc<int16_t>;
  template class Malloc<uint32_t>;
  template class Malloc<int32_t>;
  template class Malloc<uint64_t>;
  template class Malloc<int64_t>;
  template class Malloc<float>;
  template class Malloc<double>;

  template class Malloc<void *>;
  template class Malloc<char *>;
  template class Malloc<uint16_t *>;
  template class Malloc<int16_t *>;
  template class Malloc<uint32_t *>;
  template class Malloc<int32_t *>;
  template class Malloc<uint64_t *>;
  template class Malloc<int64_t *>;
  template class Malloc<float *>;
  template class Malloc<double *>;

}  // namespace Utils

#endif

//
// eof: Malloc.cc
//
