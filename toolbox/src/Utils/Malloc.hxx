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

///
/// file: Malloc.hxx
///

#pragma once

#ifndef MALLOC_dot_HXX
#define MALLOC_dot_HXX

/*\
:|:    ____            _       _             __
:|:   / ___| _     _  (_)_ __ | |_ ___ _ __ / _| __ _  ___ ___
:|:  | |   _| |_ _| |_| | '_ \| __/ _ \ '__| |_ / _` |/ __/ _ \
:|:  | |__|_   _|_   _| | | | | ||  __/ |  |  _| (_| | (_|  __/
:|:   \____||_|   |_| |_|_| |_|\__\___|_|  |_|  \__,_|\___\___|
\*/

namespace Utils {

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  using std::int64_t;
  using std::string;
  using std::mutex;
  #endif


  extern mutex MallocMutex;

  extern int64_t CountAlloc;
  extern int64_t CountFreed;
  extern int64_t AllocatedBytes;
  extern int64_t MaximumAllocatedBytes;
  extern bool    MallocDebug;

  string outBytes( size_t nb );

  /*\
  :|:   __  __       _ _
  :|:  |  \/  | __ _| | | ___   ___
  :|:  | |\/| |/ _` | | |/ _ \ / __|
  :|:  | |  | | (_| | | | (_) | (__
  :|:  |_|  |_|\__,_|_|_|\___/ \___|
  \*/

  //!
  //! Class manafin memory allocation.
  //!
  template <typename T>
  class Malloc {
  public:
    typedef T valueType;

  private:

    string      m_name;
    size_t      m_numTotValues;
    size_t      m_numTotReserved;
    size_t      m_numAllocated;
    valueType * m_pMalloc;

    Malloc(Malloc<T> const &) = delete; // blocco costruttore di copia
    Malloc<T> const & operator = (Malloc<T> &) const = delete; // blocco copia

    void allocate_internal( size_t n );
    void memory_exausted( size_t sz );

  public:

    //!
    //! Malloc object constructor
    //!
    explicit
    Malloc( string const & name );

    //!
    //! Malloc object destructor.
    //!
    ~Malloc() { hard_free(); }

    //!
    //! Allocate memory for `n` objects,
    //! raise an error if memory already allocated.
    //!
    void allocate( size_t n );

    //!
    //! Allocate memory for `n` objects,
    //! no matter if already allocated.
    //!
    void reallocate( size_t n );

    //!
    //! Free memory without deallocating pointer.
    //!
    void free(void) { m_numTotValues = m_numAllocated = 0; }

    //!
    //! Free memory deallocating pointer.
    //!
    void hard_free(void);

    //!
    //! Number of objects allocated.
    //!
    size_t size(void) const { return m_numTotValues; }

    //!
    //! Get pointer of allocated memory for `sz` objets.
    //!
    T * operator () ( size_t sz ) {
      size_t offs = m_numAllocated;
      m_numAllocated += sz;
      if ( m_numAllocated > m_numTotValues ) memory_exausted( sz );
      return m_pMalloc + offs;
    }

    T * malloc( size_t n );
    T * realloc( size_t n );

    //!
    //! `true` if you cannot get more memory pointers.
    //!
    bool is_empty() const { return m_numAllocated >= m_numTotValues; }

    //!
    //! return an error if memory is not completely used.
    //!
    void must_be_empty( char const * const where ) const;
  };

  extern template class Malloc<char>;
  extern template class Malloc<uint16_t>;
  extern template class Malloc<int16_t>;
  extern template class Malloc<uint32_t>;
  extern template class Malloc<int32_t>;
  extern template class Malloc<uint64_t>;
  extern template class Malloc<int64_t>;
  extern template class Malloc<float>;
  extern template class Malloc<double>;

  extern template class Malloc<void*>;
  extern template class Malloc<char*>;
  extern template class Malloc<uint16_t*>;
  extern template class Malloc<int16_t*>;
  extern template class Malloc<uint32_t*>;
  extern template class Malloc<int32_t*>;
  extern template class Malloc<uint64_t*>;
  extern template class Malloc<int64_t*>;
  extern template class Malloc<float*>;
  extern template class Malloc<double*>;

}

#endif

///
/// eof: Malloc.hxx
///
