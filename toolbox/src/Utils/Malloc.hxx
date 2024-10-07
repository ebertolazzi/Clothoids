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

//
// file: Malloc.hxx
//

/*\
:|:    ____            _       _             __
:|:   / ___| _     _  (_)_ __ | |_ ___ _ __ / _| __ _  ___ ___
:|:  | |   _| |_ _| |_| | '_ \| __/ _ \ '__| |_ / _` |/ __/ _ \
:|:  | |__|_   _|_   _| | | | | ||  __/ |  |  _| (_| | (_|  __/
:|:   \____||_|   |_| |_|_| |_|\__\___|_|  |_|  \__,_|\___\___|
\*/

namespace Utils {

  /*!
   * \addtogroup Malloc
   * @{
   */

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  using std::int64_t;
  using std::string;
  using std::mutex;
  #endif

  //! Global mutex for thread-safe memory operations.
  extern std::mutex MallocMutex;

  //! Global variables for tracking memory allocation statistics.
  extern int64_t CountAlloc;
  extern int64_t CountFreed;
  extern int64_t AllocatedBytes;
  extern int64_t MaximumAllocatedBytes;
  extern bool    MallocDebug;

  //! Utility function to convert byte size into a human-readable format.
  /*!
   * \param nb The number of bytes.
   * \return A string representing the size in human-readable format (KB, MB, etc.).
   */
  std::string out_bytes( size_t nb );

  /*\
  :|:   __  __       _ _
  :|:  |  \/  | __ _| | | ___   ___
  :|:  | |\/| |/ _` | | |/ _ \ / __|
  :|:  | |  | | (_| | | | (_) | (__
  :|:  |_|  |_|\__,_|_|_|\___/ \___|
  \*/

  //! Class for dynamic memory allocation of objects.
  /*!
   * This class provides custom memory management utilities for allocating,
   * freeing, and managing dynamic memory for objects of type `T`.
   *
   * \tparam T The type of objects to allocate.
   */
  template <typename T>
  class Malloc {
  public:
    //! Type alias for the type of objects managed by this allocator.
    using valueType = T;

  private:

    std::string m_name;                   //!< Name identifier for the allocated memory.
    std::size_t m_num_total_values{0};    //!< Total number of objects allocated.
    std::size_t m_num_total_reserved{0};  //!< Total reserved space.
    std::size_t m_num_allocated{0};       //!< Number of currently allocated objects.
    valueType * m_p_memory{nullptr};      //!< Pointer to the allocated memory.

    //! Internal method to allocate memory for a specified number of objects.
    void allocate_internal( std::size_t n );

    //! Handle memory exhaustion errors.
    void memory_exausted( std::size_t sz );

    //! Handle errors when attempting to pop more than allocated.
    void pop_exausted( std::size_t sz );

  public:

    //! Copy constructor is deleted.
    Malloc( Malloc<T> const & ) = delete;

    //! Assignment operator is deleted.
    Malloc<T> const & operator = ( Malloc<T> const & ) const = delete;

    //! Constructor.
    /*!
     * \param name A string identifier for the allocated memory block.
     */
    explicit
    Malloc( string name )
    : m_name(std::move(name))
    { }

    //! Destructor.
    /*!
     * Frees the allocated memory.
     */
    ~Malloc() { hard_free(); }

    //! Allocate memory for `n` objects, error if already allocated.
    /*!
     * \param n Number of objects to allocate.
     */
    void allocate( std::size_t n );

    //! Reallocate memory for `n` objects, even if already allocated.
    /*!
     * \param n Number of objects to reallocate.
     */
    void reallocate( std::size_t n );

    //! Free memory without deallocating the pointer.
    void free() { m_num_total_values = m_num_allocated = 0; }

    //! Free memory and deallocate the pointer.
    void hard_free();

    //! Get the number of allocated objects.
    /*!
     * \return Number of currently allocated objects.
     */
    size_t size() const { return m_num_total_values; }

    //! Allocate memory for `sz` objects and return the pointer.
    /*!
     * \param sz Number of objects to allocate.
     * \return Pointer to the allocated memory.
     */
    T * operator () ( std::size_t sz ) {
      size_t offs = m_num_allocated;
      m_num_allocated += sz;
      if ( m_num_allocated > m_num_total_values ) memory_exausted( sz );
      return m_p_memory + offs;
    }

    //! Free memory for `sz` objects.
    /*!
     * \param sz Number of objects to free.
     */
    void
    pop( std::size_t sz ) {
      if ( sz > m_num_allocated ) pop_exausted( sz );
      m_num_allocated -= sz;
    }

    //! Allocate memory for `n` objects.
    /*!
     * \param n Number of objects to allocate.
     * \return Pointer to the allocated memory.
     */
    T * malloc( std::size_t n );

    //! Reallocate memory for `n` objects.
    /*!
     * \param n Number of objects to reallocate.
     * \return Pointer to the reallocated memory.
     */
    T * realloc( std::size_t n );

    //! Check if the memory is fully allocated.
    /*!
     * \return `true` if all memory is allocated, `false` otherwise.
     */
    bool is_empty() const { return m_num_allocated >= m_num_total_values; }

    //! Ensure that memory is fully used.
    /*!
     * \param where Identifier for where the check is performed.
     */
    void must_be_empty( char const * where ) const;

    //! Get memory allocation information.
    /*!
     * \param where Identifier for where the information is retrieved.
     * \return A string containing information about memory allocation.
     */
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

  //! Class for fixed-size memory allocation of objects.
  /*!
   * This class manages memory for a fixed number of objects of type `T`.
   *
   * \tparam T The type of objects to allocate.
   * \tparam mem_size The fixed size of memory to allocate.
   */
  template <typename T, std::size_t mem_size>
  class MallocFixed {
  public:
    //! Type alias for the type of objects managed by this allocator.
    using valueType = T;

  private:

    std::string m_name;             //!< Name identifier for the allocated memory.
    std::size_t m_num_allocated{0}; //!< Number of currently allocated objects.
    valueType   m_data[mem_size];   //!< Array to store objects of type `T`.

  public:

    //! Copy constructor is deleted.
    MallocFixed(MallocFixed<T,mem_size> const &) = delete; // blocco costruttore di copia

    //! Assignment operator is deleted.
    MallocFixed<T,mem_size> const & operator = (MallocFixed<T,mem_size> const &) const = delete; // blocco copia

    //! Constructor.
    /*!
     * \param name A string identifier for the allocated memory block.
     */
    explicit
    MallocFixed( std::string name )
    : m_name(std::move(name))
    {}

    //! Destructor.
    ~MallocFixed() = default;

    //! Free memory without deallocating the pointer.
    void free() { m_num_allocated = 0; }

    //! Get the number of allocated objects.
    /*!
     * \return Number of objects that can be allocated.
     */
    size_t size() const { return mem_size; }

    //! Allocate memory for `sz` objects and return the pointer.
    /*!
     * \param sz Number of objects to allocate.
     * \return Pointer to the allocated memory.
     */
    T * operator () ( std::size_t sz );

    //! Free memory for `sz` objects.
    /*!
     * \param sz Number of objects to free.
     */
    void pop( std::size_t sz );

    //! Check if the memory is fully allocated.
    /*!
     * \return `true` if all memory is allocated, `false` otherwise.
     */
    bool is_empty() const { return m_num_allocated >= mem_size; }

  };

  /*! @} */

}

//
// eof: Malloc.hxx
//
