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

  extern std::mutex MallocMutex;

  extern int64_t CountAlloc;
  extern int64_t CountFreed;
  extern int64_t AllocatedBytes;
  extern int64_t MaximumAllocatedBytes;
  extern bool    MallocDebug;

  std::string out_bytes( size_t nb );

  /*\
  :|:   __  __       _ _
  :|:  |  \/  | __ _| | | ___   ___
  :|:  | |\/| |/ _` | | |/ _ \ / __|
  :|:  | |  | | (_| | | | (_) | (__
  :|:  |_|  |_|\__,_|_|_|\___/ \___|
  \*/

  //!
  //! Class for memory allocation.
  //!
  template <typename T>
  class Malloc {
  public:
    typedef T valueType;

  private:

    std::string m_name;
    std::size_t m_num_total_values;
    std::size_t m_num_total_reserved;
    std::size_t m_num_allocated;
    valueType * m_p_memory;

    void allocate_internal( std::size_t n );
    void memory_exausted( std::size_t sz );
    void pop_exausted( std::size_t sz );

  public:

    Malloc( Malloc<T> const & ) = delete; // blocco costruttore di copia
    Malloc<T> const & operator = ( Malloc<T> const & ) const = delete; // blocco copia

    //!
    //! Malloc object constructor
    //!
    explicit
    Malloc( string name )
    : m_name(std::move(name))
    , m_num_total_values(0)
    , m_num_total_reserved(0)
    , m_num_allocated(0)
    , m_p_memory(nullptr)
    { }

    //!
    //! Malloc object destructor.
    //!
    ~Malloc() { hard_free(); }

    //!
    //! Allocate memory for `n` objects,
    //! raise an error if memory already allocated.
    //!
    void allocate( std::size_t n );

    //!
    //! Allocate memory for `n` objects,
    //! no matter if already allocated.
    //!
    void reallocate( std::size_t n );

    //!
    //! Free memory without deallocating pointer.
    //!
    void free() { m_num_total_values = m_num_allocated = 0; }

    //!
    //! Free memory deallocating pointer.
    //!
    void hard_free();

    //!
    //! Number of objects allocated.
    //!
    size_t size() const { return m_num_total_values; }

    //!
    //! Get pointer of allocated memory for `sz` objets.
    //!
    T * operator () ( std::size_t sz ) {
      size_t offs = m_num_allocated;
      m_num_allocated += sz;
      if ( m_num_allocated > m_num_total_values ) memory_exausted( sz );
      return m_p_memory + offs;
    }

    //!
    //! Free pointer of allocated memory for `sz` objets.
    //!
    void
    pop( std::size_t sz ) {
      if ( sz > m_num_allocated ) pop_exausted( sz );
      m_num_allocated -= sz;
    }

    T * malloc( std::size_t n );
    T * realloc( std::size_t n );

    //!
    //! `true` if you cannot get more memory pointers.
    //!
    bool is_empty() const { return m_num_allocated >= m_num_total_values; }

    //!
    //! return an error if memory is not completely used.
    //!
    void must_be_empty( char const * where ) const;

    //!
    //! return information of memory allocations.
    //!
    std::string info( char const * where ) const;

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

  /*\
  :|:   __  __       _ _            _____ _              _
  :|:  |  \/  | __ _| | | ___   ___|  ___(_)_  _____  __| |
  :|:  | |\/| |/ _` | | |/ _ \ / __| |_  | \ \/ / _ \/ _` |
  :|:  | |  | | (_| | | | (_) | (__|  _| | |>  <  __/ (_| |
  :|:  |_|  |_|\__,_|_|_|\___/ \___|_|   |_/_/\_\___|\__,_|
  \*/

  //!
  //! Class for memory allocation.
  //!
  template <typename T, std::size_t mem_size>
  class MallocFixed {
  public:
    typedef T valueType;

  private:

    std::string m_name;
    std::size_t m_num_allocated;
    valueType   m_data[mem_size];

  public:

    MallocFixed(MallocFixed<T,mem_size> const &) = delete; // blocco costruttore di copia
    MallocFixed<T,mem_size> const & operator = (MallocFixed<T,mem_size> const &) const = delete; // blocco copia

    //!
    //! Malloc object constructor
    //!
    explicit
    MallocFixed( std::string name )
    : m_name(std::move(name))
    , m_num_allocated(0)
    {}

    //!
    //! Malloc object destructor.
    //!
    ~MallocFixed() = default;

    //!
    //! Free memory without deallocating pointer.
    //!
    void free() { m_num_allocated = 0; }

    //!
    //! Number of objects allocated.
    //!
    size_t size() const { return mem_size; }

    //!
    //! Get pointer of allocated memory for `sz` objets.
    //!
    T * operator () ( std::size_t sz ) {
      std::size_t offs = m_num_allocated;
      m_num_allocated += sz;
      UTILS_ASSERT(
        m_num_allocated <= mem_size,
        "MallocFixed<{}>::operator () ({}) -- Memory EXAUSTED\n", m_name, sz
      );
      return m_data + offs;
    }

    //!
    //! Free pointer of allocated memory for `sz` objets.
    //!
    void
    pop( std::size_t sz ) {
      UTILS_ASSERT(
        sz <= m_num_allocated,
        "MallocFixed<{}>::pop({}) -- Not enough element on Stack\n", m_name, sz
      );
      m_num_allocated -= sz;
    }

    //!
    //! `true` if you cannot get more memory pointers.
    //!
    bool is_empty() const { return m_num_allocated >= mem_size; }

  };

}

///
/// eof: Malloc.hxx
///
