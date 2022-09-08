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
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#if defined(__llvm__) || defined(__clang__)
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#endif

#include "Utils.hh"
#include <iostream>

namespace Utils {

  using std::string;
  using std::mutex;
  using std::lock_guard;
  using std::exception;
  using std::exit;
  using std::cerr;

  std::mutex MallocMutex;

  int64_t CountAlloc            = 0;
  int64_t CountFreed            = 0;
  int64_t AllocatedBytes        = 0;
  int64_t MaximumAllocatedBytes = 0;
  bool    MallocDebug           = false;

  string
  out_bytes( size_t nb ) {
    size_t Kb = nb>>10;
    size_t Mb = Kb>>10;
    size_t Gb = Mb>>10;
    if ( Gb > 0 ) {
      size_t mb = (100*(Mb & 0x3FF))/1024;
      return fmt::format( "{}Gb (+{}Mb)", Gb, mb );
    } else if ( Mb > 0 ) {
      size_t kb = (100*(Kb & 0x3FF))/1024;
      return fmt::format( "{}Mb (+{}Kb)", Mb, kb );
    } else if ( Kb > 0 ) {
      size_t b = (100*(nb & 0x3FF))/1024;
      return fmt::format( "{}Kb (+{}bytes)", Kb, b );
    }
    return fmt::format( "{} bytes", nb );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  Malloc<T>::allocate_internal( size_t n ) {
    try {
      size_t nb;
      {
        lock_guard<mutex> lock(Utils::MallocMutex);
        nb = m_num_total_reserved*sizeof(T);
        ++CountFreed; AllocatedBytes -= nb;
      }

      delete [] m_p_memory;
      m_num_total_values   = n;
      m_num_total_reserved = n + (n>>3); // 12% more values
      m_p_memory           = new T[m_num_total_reserved];

      {
        lock_guard<mutex> lock(Utils::MallocMutex);
        ++CountAlloc;
        nb = m_num_total_reserved*sizeof(T);
        AllocatedBytes += nb;
        if ( MaximumAllocatedBytes < AllocatedBytes )
          MaximumAllocatedBytes = AllocatedBytes;
      }

      if ( MallocDebug )
        fmt::print( "Allocating {} for {}\n", out_bytes( nb ), m_name );
    }
    catch ( exception const & exc ) {
      string reason = fmt::format(
        "Memory allocation failed: {}\nTry to allocate {} bytes for {}\n",
        exc.what(), n, m_name
      );
      print_trace( __LINE__, __FILE__, reason, cerr );
      exit(0);
    }
    catch (...) {
      string reason = fmt::format(
        "Memory allocation failed for {}: memory exausted\n", m_name
      );
      print_trace( __LINE__, __FILE__, reason, cerr );
      exit(0);
    }
    m_num_total_values = n;
    m_num_allocated    = 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  Malloc<T>::allocate( size_t n ) {
    UTILS_ASSERT(
      m_num_allocated == 0,
      "Malloc[{}]::allocate( {} ), try to allocate already allocated memory!\n",
      m_name, n
    );
    if ( n > m_num_total_reserved ) allocate_internal( n );
    m_num_total_values = n;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  Malloc<T>::reallocate( size_t n ) {
    if ( n > m_num_total_reserved ) allocate_internal( n );
    m_num_total_values = n;
    m_num_allocated    = 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  T *
  Malloc<T>::malloc( size_t n ) {
    UTILS_ASSERT(
      m_num_allocated == 0,
      "Malloc[{}]::malloc( {} ), try to allocate already allocated memory!\n",
      m_name, n
    );
    if ( n > m_num_total_reserved ) allocate_internal( n );
    m_num_total_values = n;
    m_num_allocated    = n;
    return m_p_memory;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  T *
  Malloc<T>::realloc( size_t n ) {
    if ( n > m_num_total_reserved ) allocate_internal( n );
    m_num_total_values = n;
    m_num_allocated    = n;
    return m_p_memory;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  Malloc<T>::hard_free() {
    if ( m_p_memory != nullptr ) {
      size_t nb;
      {
        lock_guard<mutex> lock(Utils::MallocMutex);
        nb = m_num_total_reserved*sizeof(T);
        ++CountFreed; AllocatedBytes -= nb;
      }

      if ( MallocDebug )
        fmt::print( "Freeing {} for {}\n", out_bytes( nb ), m_name );

      delete [] m_p_memory; m_p_memory = nullptr;
      m_num_total_values   = 0;
      m_num_total_reserved = 0;
      m_num_allocated      = 0;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  Malloc<T>::memory_exausted( size_t sz ) {
    string reason = fmt::format(
      "Malloc<{}>::operator () ({}) -- Memory EXAUSTED\n", m_name, sz
    );
    print_trace( __LINE__, __FILE__, reason, cerr );
    exit(0);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  Malloc<T>::pop_exausted( size_t sz ) {
    string reason = fmt::format(
      "Malloc<{}>::pop({}) -- Not enough element on Stack\n", m_name, sz
    );
    print_trace( __LINE__, __FILE__, reason, cerr );
    exit(0);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  Malloc<T>::must_be_empty( char const * where ) const {
    if ( m_num_allocated < m_num_total_values ) {
      string tmp = fmt::format(
        "in {} {}: not fully used!\nUnused: {} values\n",
        m_name, where, m_num_total_values - m_num_allocated
      );
      print_trace( __LINE__,__FILE__, tmp, cerr );
    }
    if ( m_num_allocated > m_num_total_values ) {
      string tmp = fmt::format(
        "in {} {}: too much used!\nMore used: {} values\n",
        m_name, where, m_num_allocated - m_num_total_values
      );
      print_trace( __LINE__,__FILE__, tmp, cerr );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  std::string
  Malloc<T>::info( char const * where ) const {
    std::size_t diff = m_num_allocated > m_num_total_values ?
                       m_num_allocated - m_num_total_values :
                       m_num_total_values - m_num_allocated;
    return fmt::format(
      "in {} {}\n"
      "Allocated:  {}\n"
      "Reserved:   {}\n"
      "Difference: {} [|A-R|]\n",
      m_name, where, m_num_allocated, m_num_total_values, diff
    );
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

  template class Malloc<void*>;
  template class Malloc<char*>;
  template class Malloc<uint16_t*>;
  template class Malloc<int16_t*>;
  template class Malloc<uint32_t*>;
  template class Malloc<int32_t*>;
  template class Malloc<uint64_t*>;
  template class Malloc<int64_t*>;
  template class Malloc<float*>;
  template class Malloc<double*>;

} // end namespace lapack_wrapper

#endif

///
/// eof: Malloc.cc
///
