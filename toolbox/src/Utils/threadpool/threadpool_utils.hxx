namespace threadpool {

  template <class Iterator>
  struct is_forward_iterator {
    static const bool value = std::is_base_of<
      std::forward_iterator_tag,
      typename std::iterator_traits<Iterator>::iterator_category
    >::value;
  };

  template <class Iterator>
  struct is_input_iterator {
    static const bool value = std::is_base_of<
      std::input_iterator_tag,
      typename std::iterator_traits<Iterator>::iterator_category
    >::value;
  };

  #ifndef __ICC
    #define THREADPOOL_IMPL_EXPRESSION_CHECKER(name, expression)           \
    template<class T>                                                      \
    class name {                                                           \
      /** @cond Doxygen_Suppress */                                        \
      template<class U>                                                    \
      static auto test(int) -> decltype(expression, std::true_type());     \
      template<class>                                                      \
      static std::false_type test(...);                                    \
      /** @endcond */                                                      \
    public:                                                                \
      static const bool value =                                            \
        decltype(test<typename std::remove_reference<T>::type>(0))::value; \
    };
  #else // Work around Intel compiler bug
    #define THREADPOOL_IMPL_EXPRESSION_CHECKER(name, expression)           \
    template<class T>                                                      \
    class name {                                                           \
      template<class U>                                                    \
      static constexpr decltype(expression, 1) test(bool) { return 1; };   \
      template<class>                                                      \
      static constexpr int test(...) { return 0; }                         \
    public:                                                                \
      enum { value = test<typename std::remove_reference<T>::type>(0) };   \
    };
  #endif

  THREADPOOL_IMPL_EXPRESSION_CHECKER(is_callable, std::declval<U>()());

  /**
   * An iterator for integral values
   */
  template <class I>
  class IntegralIterator : public std::iterator<std::random_access_iterator_tag, I, I, I, I> {
    using Base = std::iterator<std::random_access_iterator_tag, I, I, I, I>;
    I m_i;
  public:
    IntegralIterator() : m_i(0) { }
    IntegralIterator(I i) : m_i(i) { }
    IntegralIterator( IntegralIterator const & i ) : m_i(i.i) { }
    typename Base::value_type operator*() const { return m_i; }
    IntegralIterator& operator++() { ++m_i; return *this; }
    IntegralIterator operator++(int) { return IntegralIterator(m_i++); }
    template<class T> IntegralIterator& operator += ( T const & t ) { m_i += t; return *this; }
    template<class T> IntegralIterator& operator -= ( T const & t ) { m_i += t; return *this; }
    IntegralIterator operator + ( I const & t ) const { return IntegralIterator(m_i + t); }
    IntegralIterator operator - ( I const & t ) const { return IntegralIterator(m_i - t); }
    typename Base::difference_type operator-(const IntegralIterator& t) const { return m_i - t.i; }
    bool operator == ( IntegralIterator const & t ) const { return m_i == t.i; }
    bool operator != ( IntegralIterator const & t ) const { return m_i != t.i; }
    bool operator <  ( IntegralIterator const & t ) const { return m_i  < t.i; }
    bool operator <= ( IntegralIterator const & t ) const { return m_i <= t.i; }
    bool operator >= ( IntegralIterator const & t ) const { return m_i >= t.i; }
    bool operator >  ( IntegralIterator const & t ) const { return m_i  > t.i; }
  };

  /**
   * Dereference a pointer or iterator.
   *
   * In the thread pool, we have to:
   *
   * 1) lock access to the iterator so other threads do not interfere;
   * 2) derefence the iterator, taking a value from it;
   * 3) increment the iterator;
   * 4) unlock access to the iterator so other threads can access it.
   * 5) process the values taken from the queue.
   *
   * For input iterators, the standard says that "*a is convertible
   * to T", with a an iterator and T the value type. For the
   * increment operator ++r (with iterator r) it says "post: any
   * copies of the previous value of r are no longer required to be
   * dereferenceable". This means we can not copy the iterator in
   * step 2), since the increment in step 3) or later on the
   * increment by other threads would invalidate it. Also, we can
   * not derefernce the iterator and take the address of the
   * result. For an istream iterator, the value is overwritten when
   * the iterator is incremented. The only solution is to
   * dereference the iterator, take its value type and store a copy
   * of it before we increment the iterator.
   *
   * This is fine for primitive value types, but bad for class types
   * where copying is expensive. If we have to create a copy of each
   * container element before processing it, the speed advantage of
   * using multiple threads might be used only for all this copying.
   *
   * In the case of array pointers and container elements, copying
   * them should be unnecessary, because we can just take the
   * address of the element in step 2) above before we increment the
   * iterator, and the address of the element will remain valid when
   * the iterator is incremented.
   *
   * We handle two cases separately:
   *
   * Array pointers and container iterators are forward
   * iterators. For forward iterators, the standard specifies that
   * they offer the multi-pass guarantee, which means we can copy
   * the iterators and use both copies independently to access the
   * values.  We use the iterator traits to find out whether the
   * iterators are forward iterators, and if they are, in pass 2)
   * above we do not dereference the iterator but instead we store a
   * copy of the iterator itself and dereference it at the beginning
   * of step 5).
   *
   * For all other iterator types, just make a copy in step 2). Take
   * care that for move iterators, the copy is a moving copy.
   */
  template <class X, class Enable = void>
  struct iterval_traits;

  /**
   *  Forward iterators (e.g. all container iterators)
   *
   *  - Store a copy of the iterator.
   *  - Copying is a no-op, pass through iterator value.
   *  - Pass to function by dereferencing. Note that for move iterators,
   *    dereferencing yields a rvalue reference, and the dereferenced
   *    value may be moved by the receiving function, which is the right
   *    behaviour for move iterators.
   */
  template <class X>
  struct iterval_traits<X,typename std::enable_if<is_forward_iterator<X>::value>::type>
  {
    using type = X;
    static X const & copy( X const & r ) { return r; }
    static auto pass( X && v ) -> decltype(*std::move(v)) { return *std::move(v); }
    //static const char* name() { return "forward"; }
  };

  /**
   * Generic input iterators (e.g. input stream iterators)
   *
   * - Store a copy of the dereferenced value of the iterator.
   * - Copying involves derefencing and making a copy. Note that for
   *   move iterators, dereferencing yields an rvalue reference, and the
   *   dereferenced value may be moved, which is the right behaviour for
   *   move iterators.
   * - Pass to the function by moving the stored value.
   */
  template <class X>
  struct iterval_traits<X,typename std::enable_if<!is_forward_iterator<X>::value && is_input_iterator<X>::value>::type>
  {
    using type = typename std::decay<decltype(*std::declval<X>())>::type;
    static auto copy( X const & r ) -> decltype(*r) { return *r; } // This may move if *r is an rvalue reference
    static type&& pass( type && v ) { return std::move(v); }
    //static const char* name() { return "input"; }
  };

  /**
   * Generic output operators (e.g. output stream iterators,back inserters)
   *
   * No use as source of algorithms.
   */
  template <class X>
  struct iterval_traits<X,typename std::enable_if<!is_input_iterator<X>::value>::type>
  {
    //static const char* name() { return "output"; }
  };

  /**
   * Wrap an arbitrary thing (maybe a reference)
   */
  template <class T>
  struct Wrap {
    T value;
    template<class V> Wrap(V&& _value) : value(std::forward<V>(_value)) { }
  };

  /**
   * Call some function at the end of the current block
   */
  template <class Destructor>
  class at_scope_exit_impl {
    Destructor m_destructor;
    bool       m_active;
    at_scope_exit_impl( at_scope_exit_impl const & ) = delete;
    at_scope_exit_impl & operator=( at_scope_exit_impl const & ) = delete;
  public:
    at_scope_exit_impl() : m_active(false) { }

    explicit
    at_scope_exit_impl( Destructor&& destructor )
    : m_destructor(std::forward<Destructor>(destructor))
    , m_active(true)
    { }

    explicit
    at_scope_exit_impl( Destructor const & destructor )
    : m_destructor(destructor)
    , m_active(true)
    { }

    at_scope_exit_impl(at_scope_exit_impl&& x)
    : m_destructor(std::move(x.m_destructor))
    , m_active(x.m_active)
    { x.m_active = false; }

    at_scope_exit_impl&
    operator=(at_scope_exit_impl&& x) {
      m_destructor = std::move(x.m_destructor);
      m_active     = x.m_active;
      x.m_active   = false;
    }

    ~at_scope_exit_impl() { if (m_active) m_destructor(); }
  };

  /**
   * Create a variable that when destructed at the end of the scope
   * executes a destructor function.
   *
   * \tparam Destructor&& destructor
   *         The destructor function, maybe a lambda function.
   *
   * Use like this:
   *
   * static int a = 0;
   *
   * { // Enter scope
   *     ++a;
   *     auto x1 = at_scope_exit([&](){ --a; }
   *     // Do something, possibly throwing an exception
   * } // x1 goes out of scope, 'delete a' is called.
   */
  template<class Destructor>
  auto at_scope_exit(Destructor&& destructor)
  -> at_scope_exit_impl<Destructor>
  { return at_scope_exit_impl<Destructor>(std::forward<Destructor>(destructor)); }

  template<class Destructor>
  auto at_scope_exit(Destructor const & destructor)
	-> at_scope_exit_impl<const Destructor&>
  { return at_scope_exit_impl<Destructor const &>(destructor); }

} // End of namespace threadpool
