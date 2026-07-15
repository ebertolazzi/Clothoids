/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2013                                                      |
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
// file: GenericContainer.cc
//

#ifdef _MSC_VER
#pragma warning( disable : 4661 )
#endif

#include "GenericContainer/GenericContainer.hh"
#include <iomanip>
#include <cmath>

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#endif

#if __cplusplus >= 201103L || ( defined( _MSC_VER ) && _MSC_VER >= 1900 )
#else
#error This library needs at least a C++11 compliant compiler
#endif

#include <regex>
#include <fstream>

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define CHECK_RESIZE( V, I ) \
  if ( ( V ).size() <= ( I ) ) ( V ).resize( ( I ) + 1 )

using GC_namespace::real_type;
using std::fpclassify;

// static
// inline
// bool isZero( real_type x )
//{ return FP_ZERO == fpclassify(x); }

static bool isZero0( real_type const x )
{
  int const c = fpclassify( x );
  return FP_ZERO == c || FP_SUBNORMAL == c;
}

#endif

namespace GC_namespace
{

  //!
  //! precision used in printing number
  //!

#ifndef DOXYGEN_SHOULD_SKIP_THIS

  string to_string( complex_type const & v )
  {
    ostringstream data;
    data.precision( stream_number_precision );
    data << v.real();
    if ( v.imag() > 0 ) data << '+' << v.imag() << 'i';
    if ( v.imag() < 0 ) data << '-' << -v.imag() << 'i';
    return data.str();
  }

  void string_escape( ostream_type & stream, string const & s )
  {
    stream << '"';
    for ( auto const c : s )
    {
      if ( c == '"' ) { stream << "\\\""; }
      else if ( c == '\n' ) { stream << "\\n"; }
      else if ( c == '\r' ) { stream << "\\r"; }
      else if ( c == '\t' ) { stream << "\\t"; }
      else if ( c == '\v' ) { stream << "\\v"; }
      else if ( c == '\b' ) { stream << "\\b"; }
      else if ( c == '\a' ) { stream << "\\a"; }
      else if ( c == '\\' ) { stream << "\\\\"; }
      else
        stream << c;
    }
    stream << '"';
  }

  template <typename TYPE> ostream_type & operator<<( ostream_type & s, vector<TYPE> const & v )
  {
    s << '[';
    for ( TYPE const & vi : v ) s << ' ' << vi;
    s << " ]";
    return s;
  }

  template <> ostream_type & operator<<( ostream_type & s, vec_bool_type const & v )
  {
    s << '[';
    for ( bool const vi : v ) s << ( vi ? " true" : " false" );
    s << " ]";
    return s;
  }

  template <> ostream_type & operator<<( ostream_type & s, vec_complex_type const & v )
  {
    s << '[';
    for ( complex_type const & vi : v ) s << ' ' << to_string( vi );
    s << " ]";
    return s;
  }

  template ostream_type & operator<<( ostream_type & s, vec_int_type const & v );
  template ostream_type & operator<<( ostream_type & s, vec_long_type const & v );
  template ostream_type & operator<<( ostream_type & s, vec_real_type const & v );
  // template ostream_type & operator << ( ostream_type & s, vec_complex_type const & v );

#endif

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  template <typename TYPE> ostream_type & operator<<( ostream_type & s, mat_type<TYPE> const & mat )
  {
    if ( mat.num_rows() > 0 && mat.num_cols() )
    {
      for ( std::size_t i{ 0 }; i < mat.num_rows(); ++i )
      {
        s << std::setw( 8 ) << mat( i, 0 );
        for ( std::size_t j{ 1 }; j < mat.num_cols(); ++j ) s << " " << std::setw( 8 ) << mat( i, j );
        s << '\n';
      }
    }
    else
    {
      s << mat.num_rows() << " by " << mat.num_cols() << " matrix\n";
    }
    return s;
  }

  template <> ostream_type & operator<<( ostream_type & s, mat_complex_type const & m )
  {
    if ( m.num_rows() > 0 && m.num_cols() )
    {
      for ( std::size_t i{ 0 }; i < m.num_rows(); ++i )
      {
        s << std::setw( 8 ) << m( i, 0 );
        for ( std::size_t j{ 1 }; j < m.num_cols(); ++j ) s << " " << std::setw( 12 ) << to_string( m( i, j ) );
        s << '\n';
      }
    }
    else
    {
      s << m.num_rows() << " by " << m.num_cols() << " matrix\n";
    }
    return s;
  }

  template ostream_type & operator<<( ostream_type & s, mat_type<int_type> const & m );
  template ostream_type & operator<<( ostream_type & s, mat_type<long_type> const & m );
  template ostream_type & operator<<( ostream_type & s, mat_type<real_type> const & m );
  // template ostream_type & operator << ( ostream_type & s, mat_type<complex_type> const & m );

  class Pcre_for_GC
  {
    std::regex  reCompiled;
    std::smatch reMatches;

  public:
    Pcre_for_GC() : reCompiled( R"(^\s*\d+\s*(##?)(-|=|~|_|)\s*(.*)$)" ) {}

    ~Pcre_for_GC() = default;

    int exec( string_view const str_in, string_type matches[4] )
    {
      if ( string const str( str_in ); std::regex_match( str, reMatches, reCompiled ) )
      {
        for ( size_t i{ 0 }; i < reMatches.size(); ++i ) matches[i] = reMatches[i].str();
        return static_cast<int>( reMatches.size() );
      }
      return 0;
    }
  };

  static Pcre_for_GC pcre_for_GC;

#endif


  void GenericContainer::get_keys( vec_string_type & keys ) const
  {
    keys.clear();
    if ( GC_type::MAP == get_type() )
    {
      keys.reserve( _m().size() );
      for ( const auto & [fst, snd] : _m() ) keys.emplace_back( fst );
    }
  }

  string GenericContainer::get_keys() const
  {
    string res;
    if ( GC_type::MAP == get_type() )
    {
      for ( const auto & [fst, snd] : _m() )
      {
        res += fst;
        res += ", ";
      }
      if ( !_m().empty() )
      {
        res.pop_back();
        res.pop_back();
      }
    }
    return res;
  }

  GenericContainer & GenericContainer::operator=( vec_bool_type const & a )
  {
    set_vec_bool( static_cast<std::size_t>( a.size() ) );
    std::copy( a.begin(), a.end(), _v_b().begin() );
    return *this;
  }

  GenericContainer & GenericContainer::operator=( vec_int_type const & a )
  {
    set_vec_int( static_cast<std::size_t>( a.size() ) );
    std::copy( a.begin(), a.end(), _v_i().begin() );
    return *this;
  }

  GenericContainer & GenericContainer::operator=( vec_long_type const & a )
  {
    set_vec_long( a.size() );
    std::copy( a.begin(), a.end(), _v_l().begin() );
    return *this;
  }

  GenericContainer & GenericContainer::operator=( vec_real_type const & a )
  {
    set_vec_real( a.size() );
    std::copy( a.begin(), a.end(), _v_r().begin() );
    return *this;
  }

  GenericContainer & GenericContainer::operator=( vec_complex_type const & a )
  {
    set_vec_complex( a.size() );
    std::copy( a.begin(), a.end(), _v_c().begin() );
    return *this;
  }

  GenericContainer & GenericContainer::operator=( vec_string_type const & a )
  {
    set_vec_string( a.size() );
    std::copy( a.begin(), a.end(), _v_s().begin() );
    return *this;
  }

  GenericContainer & GenericContainer::operator=( mat_int_type const & a )
  {
    set_mat_int( a );
    return *this;
  }

  GenericContainer & GenericContainer::operator=( mat_long_type const & a )
  {
    set_mat_long( a );
    return *this;
  }

  GenericContainer & GenericContainer::operator=( mat_real_type const & a )
  {
    set_mat_real( a );
    return *this;
  }

  GenericContainer & GenericContainer::operator=( mat_complex_type const & a )
  {
    set_mat_complex( a );
    return *this;
  }

  template <typename TYPE>
  static string compare_vector( string_view const who, vector<TYPE> const & A, vector<TYPE> const & B )
  {
    ostringstream data;
    // controllo valori
    if ( A.size() != B.size() ) { data << who << " size: " << A.size() << " <> " << B.size() << '\n'; }
    else
    {
      auto it1 = A.begin();
      auto it2 = B.begin();
      int  i{ 0 };
      while ( it1 != A.end() && it2 != B.end() )
      {
        if ( *it1 != *it2 )
        {
          data << who << " at " << i << " values " << *it1 << " <> " << *it2 << '\n';
          break;
        }
        ++i;
        ++it1;
        ++it2;
      }
    }
    return data.str();
  }

  template <typename TYPE>
  static string compare_matrix( string_view const who, mat_type<TYPE> const & A, mat_type<TYPE> const & B )
  {
    ostringstream data;
    if ( A.num_rows() == B.num_rows() && A.num_cols() == B.num_cols() )
    {
      for ( std::size_t i{ 0 }; i < A.num_rows(); ++i )
      {
        for ( std::size_t j{ 0 }; j < A.num_cols(); ++j )
        {
          TYPE const & Aij = A( i, j );
          TYPE const & Bij = B( i, j );
          if ( Aij != Bij )
          {
            data << who << " at (" << i << "," << j << ") values " << Aij << " <> " << Bij << '\n';
            break;
          }
        }
      }
    }
    else
    {
      data << who << " size: " << A.num_rows() << " x " << A.num_cols() << " <> " << B.num_rows() << " x "
           << B.num_cols() << '\n';
    }
    return data.str();
  }

  string GenericContainer::compare_content( GenericContainer const & gc, string_view from ) const
  {
    ostringstream data;
    if ( get_type() != gc.get_type() )
    {
      data << from << "different type: " << to_string( get_type() ) << " <> " << to_string( gc.get_type() ) << '\n';
    }
    else
    {
      string tmp;
      switch ( get_type() )
      {
        case GC_type::NOTYPE: break;
        case GC_type::BOOL:
          if ( _b() != gc._b() ) data << from << "boolean: different\n";
          break;
        case GC_type::INTEGER:
          if ( _i() != gc._i() ) data << from << "integer: " << _i() << " <> " << gc._i() << '\n';
          break;
        case GC_type::LONG:
          if ( _l() != gc._l() ) data << from << "long: " << _l() << " <> " << gc._l() << '\n';
          break;
        case GC_type::REAL:
          if ( _r() != gc._r() ) data << from << "real: " << _r() << " <> " << gc._r() << '\n';
          break;
        case GC_type::POINTER:
          if ( _p() != gc._p() )
            data << from << "pointer: 0x" << std::hex << _p() << " <> 0x" << std::hex << gc._p() << '\n';
          break;
        case GC_type::STRING:
          if ( _s() != gc._s() ) data << from << "string: \"" << _s() << "\" <> \"" << gc._s() << "\"\n";
          break;
        case GC_type::COMPLEX:
          if ( _c() != gc._c() ) data << from << "complex: " << _c() << " <> " << gc._c() << '\n';
          break;
        case GC_type::VEC_POINTER:
          tmp = compare_vector( "vector of pointer", _v_p(), gc._v_p() );
          if ( !tmp.empty() ) data << from << tmp;
          break;
        case GC_type::VEC_BOOL:
          tmp = compare_vector( "vector of boolean", _v_b(), gc._v_b() );
          if ( !tmp.empty() ) data << from << tmp;
          break;
        case GC_type::VEC_INTEGER:
          tmp = compare_vector( "vector of integer", _v_i(), gc._v_i() );
          if ( !tmp.empty() ) data << from << tmp;
          break;
        case GC_type::VEC_LONG:
          tmp = compare_vector( "vector of long", _v_l(), gc._v_l() );
          if ( !tmp.empty() ) data << from << tmp;
          break;
        case GC_type::VEC_REAL:
          tmp = compare_vector( "vector of double", _v_r(), gc._v_r() );
          if ( !tmp.empty() ) data << from << tmp;
          break;
        case GC_type::VEC_COMPLEX:
          tmp = compare_vector( "vector of complex", _v_c(), gc._v_c() );
          if ( !tmp.empty() ) data << from << tmp;
          break;
        case GC_type::MAT_INTEGER:
          tmp = compare_matrix( "mat of integer", _m_i(), gc._m_i() );
          if ( !tmp.empty() ) data << from << tmp;
          break;
        case GC_type::MAT_LONG:
          tmp = compare_matrix( "mat of long", _m_l(), gc._m_l() );
          if ( !tmp.empty() ) data << from << tmp;
          break;
        case GC_type::MAT_REAL:
          tmp = compare_matrix( "mat of double", _m_r(), gc._m_r() );
          if ( !tmp.empty() ) data << from << tmp;
          break;
        case GC_type::MAT_COMPLEX:
          tmp = compare_matrix( "mat of complex", _m_c(), gc._m_c() );
          if ( !tmp.empty() ) data << from << tmp;
          break;
        case GC_type::VEC_STRING:
          tmp = compare_vector( "vector of string", _v_s(), gc._v_s() );
          if ( !tmp.empty() ) data << from << tmp;
          break;
        case GC_type::VECTOR:
          if ( _v().size() == gc._v().size() )
          {
            // controllo contenutp
            auto        it1 = _v().begin();
            auto        it2 = gc._v().begin();
            std::size_t i{ 0 };
            while ( it1 != _v().end() )
            {
              if ( string const res{ it1->compare_content( *it2, "> " ) }; !res.empty() )
              {
                data << from << "position: " << i << '\n' << res;
                break;
              }
              ++i;
              ++it1;
              ++it2;
            }
          }
          else
          {
            data << from << "vector of GC size do not match: " << _v().size() << " <> " << gc._v().size() << '\n';
          }
          break;
        case GC_type::MAP:
          if ( _m().size() == gc._m().size() )
          {
            // controllo le chiavi
            auto it1 = _m().begin();
            auto it2 = gc._m().begin();
            while ( it1 != _m().end() )
            {
              if ( it1->first == it2->first )
              {
                if ( string const res{ it1->second.compare_content( it2->second, "> " ) }; !res.empty() )
                {
                  data << from << "key: '" << it1->first << "'\n" << res;
                  break;
                }
              }
              else
              {
                data << from << "map of GC keys do not match: " << it1->first << " <> " << it2->first << '\n';
                break;
              }
              ++it1;
              ++it2;
            }
          }
          else
          {
            data << from << "map of GC size do not match: " << _m().size() << " <> " << gc._m().size() << '\n';
          }
          break;
      }
    }
    return data.str();
  }

  void GenericContainer::clear()
  {
    // the Box destructors release every alternative, recursing through
    // nested vectors/maps
    m_data.emplace<std::monostate>();
  }

  void GenericContainer::erase( string_view const name )
  {
    GC_assert(
      GC_type::MAP == get_type(),
      "GenericContainer::erase('{}') bad data type\n"
      "expect: {}\n"
      "but data stored is of type: {}",
      name,
      to_string( GC_type::MAP ),
      to_string( get_type() ) );
    _m().erase( string( name ) );
  }

  std::size_t GenericContainer::get_num_elements() const
  {
    switch ( get_type() )
    {
      case GC_type::POINTER:
      case GC_type::BOOL:
      case GC_type::INTEGER:
      case GC_type::LONG:
      case GC_type::REAL:
      case GC_type::COMPLEX:
      case GC_type::STRING: return 1;

      case GC_type::VEC_POINTER: return static_cast<std::size_t>( _v_p().size() );
      case GC_type::VEC_BOOL: return static_cast<std::size_t>( _v_b().size() );
      case GC_type::VEC_INTEGER: return static_cast<std::size_t>( _v_i().size() );
      case GC_type::VEC_LONG: return static_cast<std::size_t>( _v_l().size() );
      case GC_type::VEC_REAL: return static_cast<std::size_t>( _v_r().size() );
      case GC_type::VEC_COMPLEX: return static_cast<std::size_t>( _v_c().size() );
      case GC_type::VEC_STRING: return static_cast<std::size_t>( _v_s().size() );

      case GC_type::MAT_INTEGER: return static_cast<std::size_t>( _m_i().size() );
      case GC_type::MAT_LONG: return static_cast<std::size_t>( _m_l().size() );
      case GC_type::MAT_REAL: return static_cast<std::size_t>( _m_r().size() );
      case GC_type::MAT_COMPLEX: return static_cast<std::size_t>( _m_c().size() );

      case GC_type::VECTOR: return static_cast<std::size_t>( _v().size() );
      case GC_type::MAP: return static_cast<std::size_t>( _m().size() );
      case GC_type::NOTYPE: return 0;
    }
    return 0;
  }

  std::size_t GenericContainer::num_rows() const
  {
    switch ( get_type() )
    {
      case GC_type::POINTER:
      case GC_type::BOOL:
      case GC_type::INTEGER:
      case GC_type::LONG:
      case GC_type::REAL:
      case GC_type::COMPLEX:
      case GC_type::STRING:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_BOOL:
      case GC_type::VEC_INTEGER:
      case GC_type::VEC_LONG:
      case GC_type::VEC_REAL:
      case GC_type::VEC_COMPLEX:
      case GC_type::VEC_STRING:
      case GC_type::VECTOR: return 1;
      case GC_type::MAT_INTEGER: return _m_i().num_rows();
      case GC_type::MAT_LONG: return _m_l().num_rows();
      case GC_type::MAT_REAL: return _m_r().num_rows();
      case GC_type::MAT_COMPLEX: return _m_c().num_rows();
      case GC_type::MAP: return 1;
      case GC_type::NOTYPE: return 0;
    }
    return 0;
  }

  std::size_t GenericContainer::num_cols() const
  {
    switch ( get_type() )
    {
      case GC_type::POINTER:
      case GC_type::BOOL:
      case GC_type::INTEGER:
      case GC_type::LONG:
      case GC_type::REAL:
      case GC_type::COMPLEX:
      case GC_type::STRING: return 1;
      case GC_type::VEC_POINTER: return static_cast<std::size_t>( _v_p().size() );
      case GC_type::VEC_BOOL: return static_cast<std::size_t>( _v_b().size() );
      case GC_type::VEC_INTEGER: return static_cast<std::size_t>( _v_i().size() );
      case GC_type::VEC_LONG: return static_cast<std::size_t>( _v_l().size() );
      case GC_type::VEC_REAL: return static_cast<std::size_t>( _v_r().size() );
      case GC_type::VEC_COMPLEX: return static_cast<std::size_t>( _v_c().size() );
      case GC_type::VEC_STRING: return static_cast<std::size_t>( _v_s().size() );

      case GC_type::MAT_INTEGER: return _m_i().num_cols();
      case GC_type::MAT_LONG: return _m_l().num_cols();
      case GC_type::MAT_REAL: return _m_r().num_cols();
      case GC_type::MAT_COMPLEX: return _m_c().num_cols();

      case GC_type::VECTOR: return static_cast<std::size_t>( _v().size() );
      case GC_type::MAP: return static_cast<std::size_t>( _m().size() );
      case GC_type::NOTYPE: return 0;
    }
    return 0;
  }

  int GenericContainer::ck( TypeAllowed const tp ) const
  {
    if ( tp == get_type() ) return 0;       // ok
    if ( tp == GC_type::NOTYPE ) return 1;  //
    return 2;
  }

  void GenericContainer::ck( string_view const where, TypeAllowed const tp ) const
  {
    GC_assert(
      tp == get_type(),
      "{} bad data type, expect: {} but data stored is of type: {}",
      where,
      to_string( tp ),
      to_string( get_type() ) );
  }

  void GenericContainer::ck_or_set( string_view const where, TypeAllowed const tp )
  {
    if ( get_type() != GC_type::NOTYPE )
    {
      ck( where, tp );
      return;
    }
    // define the requested scalar type on an empty container (the historical
    // union left the value uninitialized here; the variant zero-initializes)
    switch ( tp )
    {
      case GC_type::POINTER: m_data.emplace<pointer_type>(); break;
      case GC_type::BOOL: m_data.emplace<bool_type>(); break;
      case GC_type::INTEGER: m_data.emplace<int_type>(); break;
      case GC_type::LONG: m_data.emplace<long_type>(); break;
      case GC_type::REAL: m_data.emplace<real_type>(); break;
      default: GC_assert( false, "{} ck_or_set() cannot default-define type: {}", where, to_string( tp ) );
    }
  }

  /*
   //      _    _ _                 _
   //     / \  | | | ___   ___ __ _| |_ ___
   //    / _ \ | | |/ _ \ / __/ _` | __/ _ \
   //   / ___ \| | | (_) | (_| (_| | ||  __/
   //  /_/   \_\_|_|\___/ \___\__,_|\__\___|
   */
  // Each allocate_* keeps the historical semantics: reuse the value when the
  // type already matches (resizing when a size is given), otherwise replace
  // the stored value with a fresh one. reset_to builds the new value before
  // touching the variant, so an allocation failure leaves the container
  // unchanged (strong guarantee -- the union version leaked its tag here).

  void GenericContainer::allocate_string()
  {
    if ( get_type() != GC_type::STRING ) reset_to<string_type>();
  }

  void GenericContainer::allocate_complex()
  {
    if ( get_type() != GC_type::COMPLEX ) reset_to<complex_type>();
  }

  void GenericContainer::allocate_vec_pointer( std::size_t const sz )
  {
    if ( get_type() != GC_type::VEC_POINTER ) reset_to<vec_pointer_type>();
    if ( sz > 0 ) _v_p().resize( sz );
  }

  GenericContainer & GenericContainer::free_pointer()
  {
    GC_assert(
      GC_type::POINTER == get_type() || GC_type::NOTYPE == get_type(),
      "GC::free_pointer() bad data type\n"
      "expect: {} but data stored is of type: {}",
      to_string( GC_type::POINTER ),
      to_string( get_type() ) );
    m_data.emplace<std::monostate>();
    return *this;
  }

  void GenericContainer::allocate_vec_bool( std::size_t const sz )
  {
    if ( get_type() != GC_type::VEC_BOOL ) reset_to<vec_bool_type>();
    if ( sz > 0 ) _v_b().resize( sz );
  }

  void GenericContainer::allocate_vec_int( std::size_t const sz )
  {
    if ( get_type() != GC_type::VEC_INTEGER ) reset_to<vec_int_type>();
    if ( sz > 0 ) _v_i().resize( sz );
  }

  void GenericContainer::allocate_vec_long( std::size_t const sz )
  {
    if ( get_type() != GC_type::VEC_LONG ) reset_to<vec_long_type>();
    if ( sz > 0 ) _v_l().resize( sz );
  }

  void GenericContainer::allocate_vec_real( std::size_t const sz )
  {
    if ( get_type() != GC_type::VEC_REAL ) reset_to<vec_real_type>();
    if ( sz > 0 ) _v_r().resize( sz );
  }

  void GenericContainer::allocate_vec_complex( std::size_t const sz )
  {
    if ( get_type() != GC_type::VEC_COMPLEX ) reset_to<vec_complex_type>();
    if ( sz > 0 ) _v_c().resize( sz );
  }

  void GenericContainer::allocate_mat_int( std::size_t const nr, std::size_t const nc )
  {
    if ( get_type() != GC_type::MAT_INTEGER )
      reset_to<mat_int_type>( nr, nc );
    else
      _m_i().resize( nr, nc );
  }

  void GenericContainer::allocate_mat_long( std::size_t const nr, std::size_t const nc )
  {
    if ( get_type() != GC_type::MAT_LONG )
      reset_to<mat_long_type>( nr, nc );
    else
      _m_l().resize( nr, nc );
  }

  void GenericContainer::allocate_mat_real( std::size_t const nr, std::size_t const nc )
  {
    if ( get_type() != GC_type::MAT_REAL )
      reset_to<mat_real_type>( nr, nc );
    else
      _m_r().resize( nr, nc );
  }

  void GenericContainer::allocate_mat_complex( std::size_t const nr, std::size_t const nc )
  {
    if ( get_type() != GC_type::MAT_COMPLEX )
      reset_to<mat_complex_type>( nr, nc );
    else
      _m_c().resize( nr, nc );
  }

  void GenericContainer::allocate_vec_string( std::size_t const sz )
  {
    if ( get_type() != GC_type::VEC_STRING ) reset_to<vec_string_type>();
    if ( sz > 0 ) _v_s().resize( sz );
  }

  void GenericContainer::allocate_vector( std::size_t const sz )
  {
    if ( get_type() != GC_type::VECTOR ) reset_to<vector_type>();
    if ( sz > 0 ) _v().resize( sz );
  }

  void GenericContainer::allocate_map()
  {
    if ( get_type() != GC_type::MAP ) reset_to<map_type>();
  }

  /*
  //   ____       _
  //  / ___|  ___| |_
  //  \___ \ / _ \ __|
  //   ___) |  __/ |_
  //  |____/ \___|\__|
  */

  pointer_type & GenericContainer::set_pointer( pointer_type const value )
  { return m_data.emplace<pointer_type>( value ); }

  bool_type & GenericContainer::set_bool( bool_type const value )
  { return m_data.emplace<bool_type>( value ); }

  int_type & GenericContainer::set_int( int_type const value )
  { return m_data.emplace<int_type>( value ); }

  long_type & GenericContainer::set_long( long_type const value )
  { return m_data.emplace<long_type>( value ); }

  real_type & GenericContainer::set_real( real_type const value )
  { return m_data.emplace<real_type>( value ); }

  complex_type & GenericContainer::set_complex( complex_type const & value )
  { return reset_to<complex_type>( value ); }

  complex_type & GenericContainer::set_complex( real_type const re, real_type const im )
  { return reset_to<complex_type>( re, im ); }

  string_type & GenericContainer::set_string( string_view const value )
  {
    allocate_string();
    return ( _s() = value );
  }

  vec_pointer_type & GenericContainer::set_vec_pointer( std::size_t const sz )
  {
    allocate_vec_pointer( sz );
    return _v_p();
  }

  vec_pointer_type & GenericContainer::set_vec_pointer( vec_pointer_type const & v )
  {
    allocate_vec_pointer( static_cast<std::size_t>( v.size() ) );
    std::copy( v.begin(), v.end(), _v_p().begin() );
    return _v_p();
  }

  vec_bool_type & GenericContainer::set_vec_bool( std::size_t const sz )
  {
    allocate_vec_bool( sz );
    return _v_b();
  }

  vec_bool_type & GenericContainer::set_vec_bool( vec_bool_type const & v )
  {
    allocate_vec_bool( static_cast<std::size_t>( v.size() ) );
    std::copy( v.begin(), v.end(), _v_b().begin() );
    return _v_b();
  }

  vec_int_type & GenericContainer::set_vec_int( std::size_t const sz )
  {
    allocate_vec_int( sz );
    return _v_i();
  }

  vec_int_type & GenericContainer::set_vec_int( vec_int_type const & v )
  {
    allocate_vec_int( static_cast<std::size_t>( v.size() ) );
    std::copy( v.begin(), v.end(), _v_i().begin() );
    return _v_i();
  }

  vec_long_type & GenericContainer::set_vec_long( std::size_t const sz )
  {
    allocate_vec_long( sz );
    return _v_l();
  }

  vec_long_type & GenericContainer::set_vec_long( vec_long_type const & v )
  {
    allocate_vec_long( static_cast<std::size_t>( v.size() ) );
    std::copy( v.begin(), v.end(), _v_l().begin() );
    return _v_l();
  }

  vec_real_type & GenericContainer::set_vec_real( std::size_t const sz )
  {
    allocate_vec_real( sz );
    return _v_r();
  }

  vec_real_type & GenericContainer::set_vec_real( vec_real_type const & v )
  {
    allocate_vec_real( static_cast<std::size_t>( v.size() ) );
    std::copy( v.begin(), v.end(), _v_r().begin() );
    return _v_r();
  }

  vec_complex_type & GenericContainer::set_vec_complex( std::size_t const sz )
  {
    allocate_vec_complex( sz );
    return _v_c();
  }

  vec_complex_type & GenericContainer::set_vec_complex( vec_complex_type const & v )
  {
    allocate_vec_complex( static_cast<std::size_t>( v.size() ) );
    std::copy( v.begin(), v.end(), _v_c().begin() );
    return _v_c();
  }

  mat_int_type & GenericContainer::set_mat_int( std::size_t const nr, std::size_t const nc )
  {
    allocate_mat_int( nr, nc );
    return _m_i();
  }

  mat_int_type & GenericContainer::set_mat_int( mat_int_type const & m )
  {
    allocate_mat_int( m.num_rows(), m.num_cols() );
    std::copy( m.begin(), m.end(), _m_i().begin() );
    return _m_i();
  }

  mat_long_type & GenericContainer::set_mat_long( std::size_t const nr, std::size_t const nc )
  {
    allocate_mat_long( nr, nc );
    return _m_l();
  }

  mat_long_type & GenericContainer::set_mat_long( mat_long_type const & m )
  {
    allocate_mat_long( m.num_rows(), m.num_cols() );
    std::copy( m.begin(), m.end(), _m_l().begin() );
    return _m_l();
  }

  mat_real_type & GenericContainer::set_mat_real( std::size_t const nr, std::size_t const nc )
  {
    allocate_mat_real( nr, nc );
    return _m_r();
  }

  mat_real_type & GenericContainer::set_mat_real( mat_real_type const & m )
  {
    allocate_mat_real( m.num_rows(), m.num_cols() );
    std::copy( m.begin(), m.end(), _m_r().begin() );
    return _m_r();
  }

  mat_complex_type & GenericContainer::set_mat_complex( std::size_t const nr, std::size_t const nc )
  {
    allocate_mat_complex( nr, nc );
    return _m_c();
  }

  mat_complex_type & GenericContainer::set_mat_complex( mat_complex_type const & m )
  {
    allocate_mat_complex( m.num_rows(), m.num_cols() );
    std::copy( m.begin(), m.end(), _m_c().begin() );
    return _m_c();
  }

  vec_string_type & GenericContainer::set_vec_string( std::size_t const sz )
  {
    allocate_vec_string( sz );
    return _v_s();
  }

  vec_string_type & GenericContainer::set_vec_string( vec_string_type const & v )
  {
    allocate_vec_string( static_cast<std::size_t>( v.size() ) );
    std::copy( v.begin(), v.end(), _v_s().begin() );
    return _v_s();
  }

  vector_type & GenericContainer::set_vector( std::size_t const sz )
  {
    allocate_vector( sz );
    return _v();
  }

  map_type & GenericContainer::set_map()
  {
    allocate_map();
    return _m();
  }

  /*
  //   ____            _
  //  |  _ \ _   _ ___| |__
  //  | |_) | | | / __| '_ \
  //  |  __/| |_| \__ \ | | |
  //  |_|    \__,_|___/_| |_|
  */
  void GenericContainer::push_bool( bool const b )
  {
    if ( get_type() == GC_type::VEC_BOOL ) { _v_b().push_back( b ); }
    else if ( get_type() == GC_type::VEC_INTEGER ) { _v_i().emplace_back( b ? 1 : 0 ); }
    else if ( get_type() == GC_type::VEC_LONG ) { _v_l().emplace_back( b ? 1 : 0 ); }
    else if ( get_type() == GC_type::VEC_REAL ) { _v_r().emplace_back( b ? 1 : 0 ); }
    else if ( get_type() == GC_type::VEC_COMPLEX ) { _v_c().emplace_back( b ? 1 : 0, 0 ); }
    else if ( get_type() == GC_type::VECTOR )
    {
      _v().resize( _v().size() + 1 );
      _v().back().set_bool( b );
    }
    else
    {
      GC_assert( false, "push_bool, bad data stored: {}", get_type_name() );
    }
  }

  void GenericContainer::push_int( int_type const i )
  {
    if ( get_type() == GC_type::VEC_INTEGER ) { _v_i().emplace_back( i ); }
    else if ( get_type() == GC_type::VEC_LONG ) { _v_l().emplace_back( static_cast<long_type>( i ) ); }
    else if ( get_type() == GC_type::VEC_REAL ) { _v_r().emplace_back( i ); }
    else if ( get_type() == GC_type::VEC_COMPLEX ) { _v_c().emplace_back( static_cast<real_type>( i ), 0 ); }
    else if ( get_type() == GC_type::VECTOR )
    {
      _v().resize( _v().size() + 1 );
      _v().back().set_int( i );
    }
    else
    {
      if ( get_type() != GC_type::VEC_INTEGER ) promote_to_vec_int();
      _v_i().emplace_back( i );
    }
  }

  void GenericContainer::push_long( long_type const l )
  {
    if ( get_type() == GC_type::VEC_LONG ) { _v_l().emplace_back( l ); }
    else if ( get_type() == GC_type::VEC_REAL ) { _v_r().emplace_back( static_cast<real_type>( l ) ); }
    else if ( get_type() == GC_type::VEC_COMPLEX ) { _v_c().emplace_back( static_cast<real_type>( l ), 0 ); }
    else if ( get_type() == GC_type::VECTOR )
    {
      _v().resize( _v().size() + 1 );
      _v().back().set_long( l );
    }
    else
    {
      if ( get_type() != GC_type::VEC_LONG ) promote_to_vec_long();
      _v_l().emplace_back( l );
    }
  }

  void GenericContainer::push_real( real_type const r )
  {
    if ( get_type() == GC_type::VEC_REAL ) { _v_r().emplace_back( r ); }
    else if ( get_type() == GC_type::VEC_COMPLEX ) { _v_c().emplace_back( r, 0 ); }
    else if ( get_type() == GC_type::VECTOR )
    {
      _v().resize( _v().size() + 1 );
      _v().back().set_real( r );
    }
    else
    {
      if ( get_type() != GC_type::VEC_REAL ) promote_to_vec_real();
      _v_r().emplace_back( r );
    }
  }

  void GenericContainer::push_complex( complex_type & c )
  {
    if ( get_type() == GC_type::VECTOR )
    {
      _v().resize( _v().size() + 1 );
      _v().back().set_complex( c );
    }
    else
    {
      if ( get_type() != GC_type::VEC_COMPLEX ) promote_to_vec_complex();
      _v_c().emplace_back( c );
    }
  }

  void GenericContainer::push_complex( real_type const re, real_type const im )
  {
    complex_type tmp( re, im );
    push_complex( tmp );
  }

  void GenericContainer::push_string( string_view const s )
  {
    if ( get_type() != GC_type::VEC_STRING ) promote_to_vector();
    if ( get_type() == GC_type::VEC_STRING ) { _v_s().emplace_back( s ); }
    else
    {
      _v().resize( _v().size() + 1 );
      _v().back().set_string( s );
    }
  }

  /*
  //    ____      _
  //   / ___| ___| |_
  //  | |  _ / _ \ __|
  //  | |_| |  __/ |_
  //   \____|\___|\__|
  */

  void * GenericContainer::get_pvoid( string_view const where ) const
  {
    ck( where, GC_type::POINTER );
    return _p();
  }

  void ** GenericContainer::get_ppvoid( string_view const where ) const
  {
    ck( where, GC_type::POINTER );
    return const_cast<void **>( &_p() );
  }

  int_type const * GenericContainer::get_int_pointer() const
  {
    switch ( get_type() )
    {
      case GC_type::INTEGER: return &_i();
      case GC_type::VEC_INTEGER: return _v_i().data();
      case GC_type::MAT_INTEGER: return _m_i().data();
      case GC_type::NOTYPE:
      case GC_type::BOOL:
      case GC_type::LONG:
      case GC_type::REAL:
      case GC_type::POINTER:
      case GC_type::STRING:
      case GC_type::COMPLEX:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_BOOL:
      case GC_type::VEC_LONG:
      case GC_type::VEC_REAL:
      case GC_type::VEC_COMPLEX:
      case GC_type::MAT_LONG:
      case GC_type::MAT_REAL:
      case GC_type::MAT_COMPLEX:
      case GC_type::VEC_STRING:
      case GC_type::VECTOR:
      case GC_type::MAP:
        GC_assert(
          false,
          "get_int_pointer, bad data type: `{}' cannot be referred as `int_type const*'",
          to_string( get_type() ) );
    }
    return nullptr;
  }

  int_type * GenericContainer::get_int_pointer()
  {
    switch ( get_type() )
    {
      case GC_type::INTEGER: return &_i();
      case GC_type::VEC_INTEGER: return _v_i().data();
      case GC_type::MAT_INTEGER: return _m_i().data();
      case GC_type::NOTYPE:
      case GC_type::BOOL:
      case GC_type::LONG:
      case GC_type::REAL:
      case GC_type::POINTER:
      case GC_type::STRING:
      case GC_type::COMPLEX:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_BOOL:
      case GC_type::VEC_LONG:
      case GC_type::VEC_REAL:
      case GC_type::VEC_COMPLEX:
      case GC_type::MAT_LONG:
      case GC_type::MAT_REAL:
      case GC_type::MAT_COMPLEX:
      case GC_type::VEC_STRING:
      case GC_type::VECTOR:
      case GC_type::MAP:
        GC_assert(
          false,
          "get_int_pointer, bad data type: `{}' cannot be referred as `int_type*'",
          to_string( get_type() ) );
    }
    return nullptr;
  }

  long_type const * GenericContainer::get_long_pointer() const
  {
    switch ( get_type() )
    {
      case GC_type::LONG: return &_l();
      case GC_type::VEC_LONG: return _v_l().data();
      case GC_type::MAT_LONG: return _m_l().data();
      case GC_type::NOTYPE:
      case GC_type::BOOL:
      case GC_type::INTEGER:
      case GC_type::REAL:
      case GC_type::POINTER:
      case GC_type::STRING:
      case GC_type::COMPLEX:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_BOOL:
      case GC_type::VEC_INTEGER:
      case GC_type::VEC_REAL:
      case GC_type::VEC_COMPLEX:
      case GC_type::MAT_INTEGER:
      case GC_type::MAT_REAL:
      case GC_type::MAT_COMPLEX:
      case GC_type::VEC_STRING:
      case GC_type::VECTOR:
      case GC_type::MAP:
        GC_assert(
          false,
          "get_long_pointer, bad data type: `{}' cannot be referred as `long_type const*'",
          to_string( get_type() ) );
    }
    return nullptr;
  }

  long_type * GenericContainer::get_long_pointer()
  {
    switch ( get_type() )
    {
      case GC_type::LONG: return &_l();
      case GC_type::VEC_LONG: return _v_l().data();
      case GC_type::MAT_LONG: return _m_l().data();
      case GC_type::NOTYPE:
      case GC_type::BOOL:
      case GC_type::INTEGER:
      case GC_type::REAL:
      case GC_type::POINTER:
      case GC_type::STRING:
      case GC_type::COMPLEX:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_BOOL:
      case GC_type::VEC_INTEGER:
      case GC_type::VEC_REAL:
      case GC_type::VEC_COMPLEX:
      case GC_type::MAT_INTEGER:
      case GC_type::MAT_REAL:
      case GC_type::MAT_COMPLEX:
      case GC_type::VEC_STRING:
      case GC_type::VECTOR:
      case GC_type::MAP:
        GC_assert(
          false,
          "get_long_pointer, bad data type: `{}' cannot be referred as `long_type*'",
          to_string( get_type() ) );
    }
    return nullptr;
  }

  real_type const * GenericContainer::get_real_pointer() const
  {
    switch ( get_type() )
    {
      case GC_type::REAL: return &_r();
      case GC_type::VEC_REAL: return _v_r().data();
      case GC_type::MAT_REAL: return _m_r().data();
      case GC_type::NOTYPE:
      case GC_type::BOOL:
      case GC_type::INTEGER:
      case GC_type::LONG:
      case GC_type::POINTER:
      case GC_type::STRING:
      case GC_type::COMPLEX:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_BOOL:
      case GC_type::VEC_INTEGER:
      case GC_type::VEC_LONG:
      case GC_type::VEC_COMPLEX:
      case GC_type::MAT_INTEGER:
      case GC_type::MAT_LONG:
      case GC_type::MAT_COMPLEX:
      case GC_type::VEC_STRING:
      case GC_type::VECTOR:
      case GC_type::MAP:
        GC_assert(
          false,
          "get_real_pointer, bad data type: `{}' cannot be referred as `real_type cont *'",
          to_string( get_type() ) );
    }
    return nullptr;
  }

  real_type * GenericContainer::get_real_pointer()
  {
    switch ( get_type() )
    {
      case GC_type::REAL: return &_r();
      case GC_type::VEC_REAL: return _v_r().data();
      case GC_type::MAT_REAL: return _m_r().data();
      case GC_type::NOTYPE:
      case GC_type::BOOL:
      case GC_type::INTEGER:
      case GC_type::LONG:
      case GC_type::POINTER:
      case GC_type::STRING:
      case GC_type::COMPLEX:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_BOOL:
      case GC_type::VEC_INTEGER:
      case GC_type::VEC_LONG:
      case GC_type::VEC_COMPLEX:
      case GC_type::MAT_INTEGER:
      case GC_type::MAT_LONG:
      case GC_type::MAT_COMPLEX:
      case GC_type::VEC_STRING:
      case GC_type::VECTOR:
      case GC_type::MAP:
        GC_assert(
          false,
          "get_real_pointer, bad data type: `{}' cannot be referred as `real_type*'",
          to_string( get_type() ) );
    }
    return nullptr;
  }

  complex_type const * GenericContainer::get_complex_pointer() const
  {
    switch ( get_type() )
    {
      case GC_type::COMPLEX: return &_c();
      case GC_type::VEC_COMPLEX: return _v_c().data();
      case GC_type::MAT_COMPLEX: return _m_c().data();
      case GC_type::NOTYPE:
      case GC_type::BOOL:
      case GC_type::INTEGER:
      case GC_type::LONG:
      case GC_type::REAL:
      case GC_type::POINTER:
      case GC_type::STRING:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_BOOL:
      case GC_type::VEC_INTEGER:
      case GC_type::VEC_LONG:
      case GC_type::VEC_REAL:
      case GC_type::MAT_INTEGER:
      case GC_type::MAT_LONG:
      case GC_type::MAT_REAL:
      case GC_type::VEC_STRING:
      case GC_type::VECTOR:
      case GC_type::MAP:
        GC_assert(
          false,
          "get_int_pointer, bad data type: `{}' cannot be referred as `complex_type const*'",
          to_string( get_type() ) );
    }
    return nullptr;
  }

  complex_type * GenericContainer::get_complex_pointer()
  {
    switch ( get_type() )
    {
      case GC_type::COMPLEX: return &_c();
      case GC_type::VEC_COMPLEX: return _v_c().data();
      case GC_type::MAT_COMPLEX: return _m_c().data();
      case GC_type::NOTYPE:
      case GC_type::BOOL:
      case GC_type::INTEGER:
      case GC_type::LONG:
      case GC_type::REAL:
      case GC_type::POINTER:
      case GC_type::STRING:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_BOOL:
      case GC_type::VEC_INTEGER:
      case GC_type::VEC_LONG:
      case GC_type::VEC_REAL:
      case GC_type::MAT_INTEGER:
      case GC_type::MAT_LONG:
      case GC_type::MAT_REAL:
      case GC_type::VEC_STRING:
      case GC_type::VECTOR:
      case GC_type::MAP:
        GC_assert(
          false,
          "get_int_pointer, bad data type: `{}' cannot be referred as `complex_type const*'",
          to_string( get_type() ) );
    }
    return nullptr;
  }

#ifndef DOXYGEN_SHOULD_SKIP_THIS

  template <> void GenericContainer::get_value( uint_type & v, string_view const where ) const
  {
    switch ( get_type() )
    {
      case GC_type::BOOL: v = _b() ? 1 : 0; break;
      case GC_type::INTEGER:
        GC_assert(
          std::in_range<uint_type>( _i() ),
          "{} in get_value(...) negative `integer` value '{}' cannot be converted into `uint_type'",
          where,
          _i() );
        v = static_cast<uint_type>( _i() );
        break;
      case GC_type::LONG:
        GC_assert(
          std::in_range<uint_type>( _l() ),
          "{} in get_value(...) negative `long` value '{}' cannot be converted into `uint_type'",
          where,
          _l() );
        v = static_cast<uint_type>( _l() );
        break;
      case GC_type::REAL:
        GC_assert(
          _r() >= 0 && GC_details::real_fits_integral<uint_type>( _r() ),
          "{} in get_value(...) negative or fractional `real` value '{}' cannot be converted into `uint_type'",
          where,
          _r() );
        v = static_cast<uint_type>( _r() );
        break;
      case GC_type::COMPLEX:
        GC_assert(
          isZero0( _c().imag() ) && GC_details::real_fits_integral<uint_type>( _c().real() ),
          "{} in get_value(...) `complex` value = {} cannot be converted into `uint_type'",
          where,
          to_string( _c() ) );
        v = static_cast<uint_type>( _c().real() );
        break;
      case GC_type::NOTYPE:
      case GC_type::POINTER:
      case GC_type::STRING:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_BOOL:
      case GC_type::VEC_INTEGER:
      case GC_type::VEC_LONG:
      case GC_type::VEC_REAL:
      case GC_type::VEC_COMPLEX:
      case GC_type::VEC_STRING:
      case GC_type::MAT_INTEGER:
      case GC_type::MAT_LONG:
      case GC_type::MAT_REAL:
      case GC_type::MAT_COMPLEX:
      case GC_type::VECTOR:
      case GC_type::MAP:
        GC_assert(
          false,
          "{}\nbad data type: `{}' cannot be converted into `uint_type'",
          where,
          to_string( get_type() ) );
    }
  }

  template <> void GenericContainer::get_value( int_type & v, string_view const where ) const
  {
    switch ( get_type() )
    {
      case GC_type::BOOL: v = _b() ? 1 : 0; break;
      case GC_type::INTEGER: v = _i(); break;
      case GC_type::LONG:
        GC_assert(
          std::in_range<int_type>( _l() ),
          "{} in get_value(...) out of range `long` value '{}' cannot be converted into `int_type'",
          where,
          _l() );
        v = static_cast<int_type>( _l() );
        break;
      case GC_type::REAL:
        GC_assert(
          GC_details::real_fits_integral<int_type>( _r() ),
          "{} in get_value(...) fractional `real` value '{}' cannot be converted into `int_type'",
          where,
          _r() );
        v = static_cast<int>( _r() );
        break;
      case GC_type::COMPLEX:
        GC_assert(
          isZero0( _c().imag() ) && GC_details::real_fits_integral<int_type>( _c().real() ),
          "{} in get_value(...) `complex` value = {} cannot be converted into `int_type'",
          where,
          to_string( _c() ) );
        v = static_cast<int>( _c().real() );
        break;
      case GC_type::NOTYPE:
      case GC_type::POINTER:
      case GC_type::STRING:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_BOOL:
      case GC_type::VEC_INTEGER:
      case GC_type::VEC_LONG:
      case GC_type::VEC_REAL:
      case GC_type::VEC_COMPLEX:
      case GC_type::VEC_STRING:
      case GC_type::MAT_INTEGER:
      case GC_type::MAT_LONG:
      case GC_type::MAT_REAL:
      case GC_type::MAT_COMPLEX:
      case GC_type::VECTOR:
      case GC_type::MAP:
        GC_assert(
          false,
          "{} in get_value(...) bad data type: `{}' cannot be converted into `int_type'",
          where,
          to_string( get_type() ) );
    }
  }

  template <> void GenericContainer::get_value( ulong_type & v, string_view const where ) const
  {
    switch ( get_type() )
    {
      case GC_type::BOOL: v = _b() ? 1 : 0; break;
      case GC_type::INTEGER:
        GC_assert(
          std::in_range<ulong_type>( _i() ),
          "{} in get_value(...) negative `integer` value '{}' cannot be converted into `ulong_type'",
          where,
          _i() );
        v = static_cast<ulong_type>( _i() );
        break;
      case GC_type::LONG:
        GC_assert(
          std::in_range<ulong_type>( _l() ),
          "{} in get_value(...) negative `long` value '{}' cannot be converted into `ulong_type'",
          where,
          _l() );
        v = static_cast<ulong_type>( _l() );
        break;
      case GC_type::REAL:
        GC_assert(
          _r() >= 0 && GC_details::real_fits_integral<ulong_type>( _r() ),
          "{} in get_value(...) negative or fractional `real` value '{}' cannot be converted into `ulong_type'",
          where,
          _r() );
        v = static_cast<ulong_type>( _r() );
        break;
      case GC_type::COMPLEX:
        GC_assert(
          isZero0( _c().imag() ) && GC_details::real_fits_integral<ulong_type>( _c().real() ),
          "{} in get_value(...) `complex` value {} cannot be converted into `ulong_type'",
          where,
          to_string( _c() ) );
        v = static_cast<ulong_type>( _c().real() );
        break;
      case GC_type::NOTYPE:
      case GC_type::POINTER:
      case GC_type::STRING:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_BOOL:
      case GC_type::VEC_INTEGER:
      case GC_type::VEC_LONG:
      case GC_type::VEC_REAL:
      case GC_type::VEC_COMPLEX:
      case GC_type::VEC_STRING:
      case GC_type::MAT_INTEGER:
      case GC_type::MAT_LONG:
      case GC_type::MAT_REAL:
      case GC_type::MAT_COMPLEX:
      case GC_type::VECTOR:
      case GC_type::MAP:
        GC_assert(
          false,
          "{} in get_value(...) bad data type: `{}' cannot be converted into `ulong_type'",
          where,
          to_string( get_type() ) );
    }
  }

  template <> void GenericContainer::get_value( long_type & v, string_view const where ) const
  {
    switch ( get_type() )
    {
      case GC_type::BOOL: v = _b() ? 1 : 0; break;
      case GC_type::INTEGER: v = static_cast<long>( _i() ); break;
      case GC_type::LONG: v = static_cast<long>( _l() ); break;
      case GC_type::REAL:
        GC_assert(
          GC_details::real_fits_integral<long_type>( _r() ),
          "{} in get_value(...) fractional `real` value '{}' cannot be converted into `long_type'",
          where,
          _r() );
        v = static_cast<long>( _r() );
        break;
      case GC_type::COMPLEX:
        GC_assert(
          isZero0( _c().imag() ) && GC_details::real_fits_integral<long_type>( _c().real() ),
          "{} in get_value(...) `complex` value = {} cannot be converted into `long_type'",
          where,
          to_string( _c() ) );
        v = static_cast<long>( _c().real() );
        break;
      case GC_type::NOTYPE:
      case GC_type::POINTER:
      case GC_type::STRING:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_BOOL:
      case GC_type::VEC_INTEGER:
      case GC_type::VEC_LONG:
      case GC_type::VEC_REAL:
      case GC_type::VEC_COMPLEX:
      case GC_type::VEC_STRING:
      case GC_type::MAT_INTEGER:
      case GC_type::MAT_LONG:
      case GC_type::MAT_REAL:
      case GC_type::MAT_COMPLEX:
      case GC_type::VECTOR:
      case GC_type::MAP:
        GC_assert(
          false,
          "{} in get_value(...) bad data type: `{}' cannot be converted into `long_type'",
          where,
          to_string( get_type() ) );
    }
  }

  template <> void GenericContainer::get_value( float & v, string_view const where ) const
  {
    switch ( get_type() )
    {
      case GC_type::BOOL: v = static_cast<float>( _b() ? 1 : 0 ); break;
      case GC_type::INTEGER: v = static_cast<float>( _i() ); break;
      case GC_type::LONG: v = static_cast<float>( _l() ); break;
      case GC_type::REAL: v = static_cast<float>( _r() ); break;
      case GC_type::COMPLEX:
        GC_assert(
          isZero0( _c().imag() ),
          "{} in get_value(...) `complex` value = {} cannot be converted into `float'",
          where,
          to_string( _c() ) );
        v = static_cast<float>( _c().real() );
        break;
      case GC_type::NOTYPE:
      case GC_type::POINTER:
      case GC_type::STRING:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_BOOL:
      case GC_type::VEC_INTEGER:
      case GC_type::VEC_LONG:
      case GC_type::VEC_REAL:
      case GC_type::VEC_COMPLEX:
      case GC_type::VEC_STRING:
      case GC_type::MAT_INTEGER:
      case GC_type::MAT_LONG:
      case GC_type::MAT_REAL:
      case GC_type::MAT_COMPLEX:
      case GC_type::VECTOR:
      case GC_type::MAP:
        GC_assert(
          false,
          "{} in get_value(...) bad data type: `{}' cannot be converted into `float'",
          where,
          to_string( get_type() ) );
    }
  }

  template <> void GenericContainer::get_value( double & v, string_view const where ) const
  {
    switch ( get_type() )
    {
      case GC_type::BOOL: v = _b() ? 1 : 0; break;
      case GC_type::INTEGER: v = static_cast<double>( _i() ); break;
      case GC_type::LONG: v = static_cast<double>( _l() ); break;
      case GC_type::REAL: v = _r(); break;
      case GC_type::COMPLEX:
        GC_assert(
          isZero0( _c().imag() ),
          "{} in get_value(...) `complex` value = {} cannot be converted into `double'",
          where,
          to_string( _c() ) );
        v = _c().real();
        break;
      case GC_type::NOTYPE:
      case GC_type::POINTER:
      case GC_type::STRING:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_BOOL:
      case GC_type::VEC_INTEGER:
      case GC_type::VEC_LONG:
      case GC_type::VEC_REAL:
      case GC_type::VEC_COMPLEX:
      case GC_type::VEC_STRING:
      case GC_type::MAT_INTEGER:
      case GC_type::MAT_LONG:
      case GC_type::MAT_REAL:
      case GC_type::MAT_COMPLEX:
      case GC_type::VECTOR:
      case GC_type::MAP:
        GC_assert(
          false,
          "{} in get_value(...) bad data type: `{}' cannot be converted into `double'",
          where,
          to_string( get_type() ) );
    }
  }

#endif

  bool GenericContainer::is_number() const
  {
    switch ( get_type() )
    {
      case GC_type::BOOL:
      case GC_type::INTEGER:
      case GC_type::LONG:
      case GC_type::REAL: return true;
      default: return false;
    }
    return false;
  }

  real_type GenericContainer::get_number( string_view const where ) const
  {
    switch ( get_type() )
    {
      case GC_type::BOOL: return _b() ? 1 : 0;
      case GC_type::INTEGER: return static_cast<real_type>( _i() );
      case GC_type::LONG: return static_cast<real_type>( _l() );
      case GC_type::REAL: return _r();
      case GC_type::NOTYPE:
      case GC_type::POINTER:
      case GC_type::COMPLEX:
      case GC_type::STRING:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_BOOL:
      case GC_type::VEC_INTEGER:
      case GC_type::VEC_LONG:
      case GC_type::VEC_REAL:
      case GC_type::VEC_COMPLEX:
      case GC_type::VEC_STRING:
      case GC_type::MAT_INTEGER:
      case GC_type::MAT_LONG:
      case GC_type::MAT_REAL:
      case GC_type::MAT_COMPLEX:
      case GC_type::VECTOR:
      case GC_type::MAP:
        GC_assert( false, "{} get_number() type: {} cannot be converted to double.", where, to_string( get_type() ) );
    }
    return 0;
  }

  complex_type GenericContainer::get_complex_number( string_view const where ) const
  {
    switch ( get_type() )
    {
      case GC_type::BOOL: return { _b() ? static_cast<real_type>( 1 ) : 0, 0 };
      case GC_type::INTEGER: return { static_cast<real_type>( _i() ), 0 };
      case GC_type::LONG: return { static_cast<real_type>( _l() ), 0 };
      case GC_type::REAL: return { _r(), 0 };
      case GC_type::COMPLEX: return _c();
      case GC_type::NOTYPE:
      case GC_type::POINTER:
      case GC_type::STRING:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_BOOL:
      case GC_type::VEC_INTEGER:
      case GC_type::VEC_LONG:
      case GC_type::VEC_REAL:
      case GC_type::VEC_COMPLEX:
      case GC_type::VEC_STRING:
      case GC_type::MAT_INTEGER:
      case GC_type::MAT_LONG:
      case GC_type::MAT_REAL:
      case GC_type::MAT_COMPLEX:
      case GC_type::VECTOR:
      case GC_type::MAP:
        GC_assert( false, "{} get_number() type: {} cannot be converted to complex.", where, to_string( get_type() ) );
    }
    return 0;
  }

  void GenericContainer::get_complex_number( real_type & re, real_type & im ) const
  {
    complex_type const tmp{ get_complex_number() };
    re = tmp.real();
    im = tmp.imag();
  }

  real_type GenericContainer::get_number_at( std::size_t const i, string_view const where ) const
  {
    switch ( get_type() )
    {
      case GC_type::VEC_BOOL: return _v_b()[i] ? 1 : 0;
      case GC_type::VEC_INTEGER: return static_cast<real_type>( _v_i()[i] );
      case GC_type::VEC_LONG: return static_cast<real_type>( _v_l()[i] );
      case GC_type::VEC_REAL: return _v_r()[i];
      case GC_type::MAT_INTEGER: return static_cast<real_type>( _m_i()[i] );
      case GC_type::MAT_LONG: return static_cast<real_type>( _m_l()[i] );
      case GC_type::MAT_REAL: return _m_r()[i];
      case GC_type::VECTOR: return _v()[i].get_number();
      case GC_type::NOTYPE:
      case GC_type::POINTER:
      case GC_type::BOOL:
      case GC_type::INTEGER:
      case GC_type::LONG:
      case GC_type::REAL:
      case GC_type::COMPLEX:
      case GC_type::STRING:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_COMPLEX:
      case GC_type::VEC_STRING:
      case GC_type::MAT_COMPLEX:
      case GC_type::MAP:
        GC_assert(
          false,
          "{} get_number_at( {} ) type: {} cannot be converted to double.",
          where,
          i,
          to_string( get_type() ) );
    }
    return 0;
  }

  complex_type GenericContainer::get_complex_number_at( std::size_t const i, string_view const where ) const
  {
    switch ( get_type() )
    {
      case GC_type::VEC_BOOL: return { static_cast<real_type>( _v_b()[i] ? 1 : 0 ), 0 };
      case GC_type::VEC_INTEGER: return { static_cast<real_type>( _v_i()[i] ), 0 };
      case GC_type::VEC_LONG: return { static_cast<real_type>( _v_l()[i] ), 0 };
      case GC_type::VEC_REAL: return { _v_r()[i], 0 };
      case GC_type::VEC_COMPLEX: return _v_c()[i];
      case GC_type::MAT_INTEGER: return { static_cast<real_type>( _m_i()[i] ), 0 };
      case GC_type::MAT_LONG: return { static_cast<real_type>( _m_l()[i] ), 0 };
      case GC_type::MAT_REAL: return { _m_r()[i], 0 };
      case GC_type::MAT_COMPLEX: return _m_c()[i];
      case GC_type::VECTOR: return _v()[i].get_complex_number();
      case GC_type::NOTYPE:
      case GC_type::POINTER:
      case GC_type::BOOL:
      case GC_type::INTEGER:
      case GC_type::LONG:
      case GC_type::REAL:
      case GC_type::COMPLEX:
      case GC_type::STRING:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_STRING:
      case GC_type::MAP:
        GC_assert(
          false,
          "{} get_complex_number_at( {} ) type: {} cannot be converted to complex.",
          where,
          i,
          to_string( get_type() ) );
    }
    return 0;
  }

  void GenericContainer::get_complex_number_at(
    std::size_t const i,
    real_type &       re,
    real_type &       im,
    string_view const where ) const
  {
    complex_type const tmp{ get_complex_number_at( i, where ) };
    re = tmp.real();
    im = tmp.imag();
  }

  bool_type & GenericContainer::get_bool( string_view const where )
  {
    ck_or_set( where, GC_type::BOOL );
    return _b();
  }

  bool_type const & GenericContainer::get_bool( string_view const where ) const
  {
    ck( where, GC_type::BOOL );
    return _b();
  }

  int_type & GenericContainer::get_int( string_view const where )
  {
    ck_or_set( where, GC_type::INTEGER );
    return _i();
  }

  int_type const & GenericContainer::get_int( string_view const where ) const
  {
    ck( where, GC_type::INTEGER );
    return _i();
  }

  long_type & GenericContainer::get_long( string_view const where )
  {
    ck_or_set( where, GC_type::LONG );
    return _l();
  }

  long_type const & GenericContainer::get_long( string_view const where ) const
  {
    ck( where, GC_type::LONG );
    return _l();
  }

  int_type GenericContainer::get_as_int( string_view const where ) const
  {
    int_type res;
    this->get_value<int_type>( res, where );
    return res;
  }

  uint_type GenericContainer::get_as_uint( string_view const where ) const
  {
    uint_type res;
    this->get_value<uint_type>( res, where );
    return res;
  }

  long_type GenericContainer::get_as_long( string_view const where ) const
  {
    long_type res;
    this->get_value<long_type>( res, where );
    return res;
  }

  ulong_type GenericContainer::get_as_ulong( string_view const where ) const
  {
    ulong_type res;
    this->get_value<ulong_type>( res, where );
    return res;
  }

  real_type & GenericContainer::get_real( string_view const where )
  {
    ck_or_set( where, GC_type::REAL );
    return _r();
  }

  real_type const & GenericContainer::get_real( string_view const where ) const
  {
    ck( where, GC_type::REAL );
    return _r();
  }

  complex_type & GenericContainer::get_complex( string_view const where )
  {
    if ( get_type() == GC_type::NOTYPE ) return reset_to<complex_type>();
    ck( where, GC_type::COMPLEX );
    return _c();
  }

  complex_type const & GenericContainer::get_complex( string_view const where ) const
  {
    ck( where, GC_type::COMPLEX );
    return _c();
  }

  string_type & GenericContainer::get_string( string_view const where )
  {
    if ( get_type() == GC_type::NOTYPE ) return reset_to<string_type>();
    ck( where, GC_type::STRING );
    return _s();
  }

  string const & GenericContainer::get_string( string_view const where ) const
  {
    ck( where, GC_type::STRING );
    return _s();
  }

  vector_type & GenericContainer::get_vector( string_view const where )
  {
    ck( where, GC_type::VECTOR );
    return _v();
  }

  vector_type const & GenericContainer::get_vector( string_view const where ) const
  {
    ck( where, GC_type::VECTOR );
    return _v();
  }

  vec_pointer_type & GenericContainer::get_vec_pointer( string_view const where )
  {
    ck( where, GC_type::VEC_POINTER );
    return _v_p();
  }

  vec_pointer_type const & GenericContainer::get_vec_pointer( string_view const where ) const
  {
    ck( where, GC_type::VEC_POINTER );
    return _v_p();
  }

  vec_bool_type & GenericContainer::get_vec_bool( string_view const where )
  {
    ck( where, GC_type::VEC_BOOL );
    return _v_b();
  }

  vec_bool_type const & GenericContainer::get_vec_bool( string_view const where ) const
  {
    ck( where, GC_type::VEC_BOOL );
    return _v_b();
  }

  vec_int_type & GenericContainer::get_vec_int( string_view const where )
  {
    if ( get_type() == GC_type::NOTYPE ) set_vec_int();
    if ( get_type() == GC_type::VEC_BOOL ) promote_to_vec_int();
    ck( where, GC_type::VEC_INTEGER );
    return _v_i();
  }

  vec_int_type const & GenericContainer::get_vec_int( string_view const where ) const
  {
    ck( where, GC_type::VEC_INTEGER );
    return _v_i();
  }

  vec_long_type & GenericContainer::get_vec_long( string_view const where )
  {
    if ( get_type() == GC_type::NOTYPE ) set_vec_long();
    if ( get_type() == GC_type::VEC_BOOL || get_type() == GC_type::VEC_INTEGER ) promote_to_vec_long();
    ck( where, GC_type::VEC_LONG );
    return _v_l();
  }

  vec_long_type const & GenericContainer::get_vec_long( string_view const where ) const
  {
    ck( where, GC_type::VEC_LONG );
    return _v_l();
  }

  vec_real_type & GenericContainer::get_vec_real( string_view const where )
  {
    if ( get_type() == GC_type::NOTYPE ) set_vec_real();
    if ( get_type() == GC_type::VEC_BOOL || get_type() == GC_type::VEC_INTEGER || get_type() == GC_type::VEC_LONG )
      promote_to_vec_real();
    ck( where, GC_type::VEC_REAL );
    return _v_r();
  }

  vec_real_type const & GenericContainer::get_vec_real( string_view const where ) const
  {
    ck( where, GC_type::VEC_REAL );
    return _v_r();
  }

  vec_complex_type & GenericContainer::get_vec_complex( string_view const where )
  {
    if ( get_type() == GC_type::NOTYPE ) set_vec_complex();
    if (
      get_type() == GC_type::VEC_BOOL || get_type() == GC_type::VEC_INTEGER || get_type() == GC_type::VEC_LONG ||
      get_type() == GC_type::VEC_REAL )
      promote_to_vec_complex();
    ck( where, GC_type::VEC_COMPLEX );
    return _v_c();
  }

  vec_complex_type const & GenericContainer::get_vec_complex( string_view const where ) const
  {
    ck( where, GC_type::VEC_COMPLEX );
    return _v_c();
  }

  mat_int_type & GenericContainer::get_mat_int( string_view const where )
  {
    if ( get_type() == GC_type::NOTYPE ) set_mat_int();
    if (
      get_type() == GC_type::VEC_BOOL || get_type() == GC_type::VEC_INTEGER || get_type() == GC_type::VEC_LONG ||
      get_type() == GC_type::VEC_REAL )
      promote_to_mat_int();
    ck( where, GC_type::MAT_INTEGER );
    return _m_i();
  }

  mat_int_type const & GenericContainer::get_mat_int( string_view const where ) const
  {
    ck( where, GC_type::MAT_INTEGER );
    return _m_i();
  }

  mat_long_type & GenericContainer::get_mat_long( string_view const where )
  {
    if ( get_type() == GC_type::NOTYPE ) set_mat_long();
    if (
      get_type() == GC_type::VEC_BOOL || get_type() == GC_type::VEC_INTEGER || get_type() == GC_type::VEC_LONG ||
      get_type() == GC_type::VEC_REAL || get_type() == GC_type::MAT_INTEGER )
      promote_to_mat_long();
    ck( where, GC_type::MAT_LONG );
    return _m_l();
  }

  mat_long_type const & GenericContainer::get_mat_long( string_view const where ) const
  {
    ck( where, GC_type::MAT_LONG );
    return _m_l();
  }

  mat_real_type & GenericContainer::get_mat_real( string_view const where )
  {
    if ( get_type() == GC_type::NOTYPE ) set_mat_real();
    if (
      get_type() == GC_type::VEC_BOOL || get_type() == GC_type::VEC_INTEGER || get_type() == GC_type::VEC_LONG ||
      get_type() == GC_type::VEC_REAL || get_type() == GC_type::MAT_INTEGER || get_type() == GC_type::MAT_LONG )
      promote_to_mat_real();
    ck( where, GC_type::MAT_REAL );
    return _m_r();
  }

  mat_real_type const & GenericContainer::get_mat_real( string_view const where ) const
  {
    ck( where, GC_type::MAT_REAL );
    return _m_r();
  }

  mat_complex_type & GenericContainer::get_mat_complex( string_view const where )
  {
    if ( get_type() == GC_type::NOTYPE ) set_mat_complex();
    if (
      get_type() == GC_type::VEC_BOOL || get_type() == GC_type::VEC_INTEGER || get_type() == GC_type::VEC_LONG ||
      get_type() == GC_type::VEC_REAL || get_type() == GC_type::MAT_REAL || get_type() == GC_type::VEC_COMPLEX )
      promote_to_mat_complex();
    ck( where, GC_type::MAT_COMPLEX );
    return _m_c();
  }

  mat_complex_type const & GenericContainer::get_mat_complex( string_view const where ) const
  {
    ck( where, GC_type::MAT_COMPLEX );
    return _m_c();
  }

  vec_string_type & GenericContainer::get_vec_string( string_view const where )
  {
    ck( where, GC_type::VEC_STRING );
    return _v_s();
  }

  vec_string_type const & GenericContainer::get_vec_string( string_view const where ) const
  {
    ck( where, GC_type::VEC_STRING );
    return _v_s();
  }

  map_type & GenericContainer::get_map( string_view const where )
  {
    ck( where, GC_type::MAP );
    return _m();
  }

  map_type const & GenericContainer::get_map( string_view const where ) const
  {
    ck( where, GC_type::MAP );
    return _m();
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------

  bool_type GenericContainer::get_map_bool( string_view const key, string_view const where ) const
  {
    GC_assert( this->exists( key ), "{} key: `{}` is missing", where, key );
    return this->_m().at( string_type( key ) ).get_bool( where );
  }

  bool_type GenericContainer::get_map_bool( std::initializer_list<string> const args ) const
  {
    string msg = "{ ";
    for ( string const & key : args )
    {
      msg += key;
      msg += ", ";
      if ( this->exists( key ) )
      {
        if ( _m().at( key ).get_type() == GC_type::BOOL ) return this->get_map_bool( key );
      }
    }
    msg.pop_back();
    msg.pop_back();
    msg += " }";
    GC_assert( false, "keys {} of type `bool` missing in the map", msg );
    return false;
  }

  bool_type GenericContainer::get_map_bool( vec_string_type const & keys, string_view const where ) const
  {
    string const who{ must_exists( keys, where ) };
    return this->_m().at( who ).get_bool( where );
  }

  int_type GenericContainer::get_map_int( string_view const key, string_view const where ) const
  {
    GC_assert( this->exists( key ), "{} key: `{}` is missing", where, key );
    return this->_m().at( string_type( key ) ).get_as_int( where );
  }

  int_type GenericContainer::get_map_int( std::initializer_list<string> const args ) const
  {
    string msg = "{ ";
    for ( string const & key : args )
    {
      msg += key;
      msg += ", ";
      if ( this->exists( key ) )
      {
        if ( _m().at( key ).get_type() == GC_type::INTEGER ) return this->get_map_int( key );
      }
    }
    msg.pop_back();
    msg.pop_back();
    msg += " }";
    GC_assert( false, "keys {} of type `int_type` missing in the map", msg );
    return false;
  }

  int_type GenericContainer::get_map_int( vec_string_type const & keys, string_view const where ) const
  {
    string const who{ must_exists( keys, where ) };
    return this->_m().at( who ).get_as_int( where );
  }

  real_type GenericContainer::get_map_number( string_view const key, string_view const where ) const
  {
    GC_assert( this->exists( key ), "{} key: `{}` is missing", where, key );
    return this->_m().at( string_type( key ) ).get_number( where );
  }

  real_type GenericContainer::get_map_number( std::initializer_list<string> const args ) const
  {
    string msg = "{ ";
    for ( string const & key : args )
    {
      msg += key;
      msg += ", ";
      if ( this->exists( key ) )
      {
        if ( _m().at( key ).is_number() ) return this->get_map_number( key );
      }
    }
    msg.pop_back();
    msg.pop_back();
    msg += " }";
    GC_assert( false, "keys {} of type `real_type` missing in the map", msg );
    return false;
  }

  real_type GenericContainer::get_map_number( vec_string_type const & keys, string_view const where ) const
  {
    string const who{ must_exists( keys, where ) };
    return this->_m().at( who ).get_number( where );
  }

  string const & GenericContainer::get_map_string( string_view const key, string_view const where ) const
  {
    GC_assert( this->exists( key ), "{} key: `{}` is missing", where, key );
    return this->_m().at( string_type( key ) ).get_string( where );
  }

  string const & GenericContainer::get_map_string( std::initializer_list<string> const args ) const
  {
    string msg = "{ ";
    for ( string const & key : args )
    {
      msg += key;
      msg += ", ";
      if ( this->exists( key ) )
      {
        if ( _m().at( key ).get_type() == GC_type::STRING ) return this->get_map_string( key );
      }
    }
    msg.pop_back();
    msg.pop_back();
    msg += " }";
    GC_assert( false, "keys {} of type `string` missing in the map", msg );
    static string empty;
    return empty;
  }

  string const & GenericContainer::get_map_string( vec_string_type const & keys, string_view const where ) const
  {
    string const who{ must_exists( keys, where ) };
    return this->_m().at( who ).get_string( where );
  }

  vec_real_type const & GenericContainer::get_map_vec_real( string_view const key, string_view const where ) const
  {
    GC_assert( this->exists( key ), "{} key: `{}` is missing", where, key );
    return this->_m().at( string_type( key ) ).get_vec_real( where );
  }

  vec_real_type const & GenericContainer::get_map_vec_real( std::initializer_list<string> const args ) const
  {
    string msg = "{ ";
    for ( string const & key : args )
    {
      msg += key;
      msg += ", ";
      if ( this->exists( key ) )
      {
        if ( _m().at( key ).get_type() == GC_type::VEC_REAL ) return this->get_map_vec_real( key );
      }
    }
    msg.pop_back();
    msg.pop_back();
    msg += " }";
    GC_assert( false, "keys {} of type `vec_real_type` missing in the map", msg );
    static vec_real_type empty;
    return empty;
  }

  vec_real_type const & GenericContainer::get_map_vec_real( vec_string_type const & keys, string_view const where )
    const
  {
    string const who{ must_exists( keys, where ) };
    return this->_m().at( who ).get_vec_real( where );
  }

  vec_complex_type const & GenericContainer::get_map_vec_complex( string_view const key, string_view const where ) const
  {
    GC_assert( this->exists( key ), "{} key: `{}` is missing", where, key );
    return this->_m().at( string_type( key ) ).get_vec_complex( where );
  }

  vec_complex_type const & GenericContainer::get_map_vec_complex( std::initializer_list<string> const args ) const
  {
    string msg = "{ ";
    for ( string const & key : args )
    {
      msg += key;
      msg += ", ";
      if ( this->exists( key ) )
      {
        if ( _m().at( key ).get_type() == GC_type::VEC_COMPLEX ) return this->get_map_vec_complex( key );
      }
    }
    msg.pop_back();
    msg.pop_back();
    msg += " }";
    GC_assert( false, "keys {} of type `vec_complex_type` missing in the map", msg );
    static vec_complex_type empty;
    return empty;
  }

  vec_complex_type const & GenericContainer::get_map_vec_complex(
    vec_string_type const & keys,
    string_view const       where ) const
  {
    string const who{ must_exists( keys, where ) };
    return this->_m().at( who ).get_vec_complex( where );
  }

  vec_string_type const & GenericContainer::get_map_vec_string( string_view const key, string_view const where ) const
  {
    GC_assert( this->exists( key ), "{} key: `{}` is missing", where, key );
    return this->_m().at( string_type( key ) ).get_vec_string( where );
  }

  vec_string_type const & GenericContainer::get_map_vec_string( std::initializer_list<string> const args ) const
  {
    string msg = "{ ";
    for ( string const & key : args )
    {
      msg += key;
      msg += ", ";
      if ( this->exists( key ) )
      {
        if ( _m().at( key ).get_type() == GC_type::VEC_STRING ) return this->get_map_vec_string( key );
      }
    }
    msg.pop_back();
    msg.pop_back();
    msg += " }";
    GC_assert( false, "keys {} missing in the map", msg );
    static vec_string_type empty;
    return empty;
  }

  vec_string_type const & GenericContainer::get_map_vec_string( vec_string_type const & keys, string_view const where )
    const
  {
    string const who{ must_exists( keys, where ) };
    return this->_m().at( who ).get_vec_string( where );
  }

  // --------------------------------------------------------------

  bool GenericContainer::exists( string_view const s ) const
  {
    if ( get_type() != GC_type::MAP ) return false;
    auto const iv{ _m().find( string_type( s ) ) };
    return iv != _m().end();
  }

  bool GenericContainer::exists( vec_string_type const & vs ) const
  {
    if ( get_type() != GC_type::MAP ) return false;
    return std::any_of(
      vs.begin(),
      vs.end(),
      [this]( string_type const & s ) { return _m().find( s ) != _m().end(); } );
  }
  string GenericContainer::must_exists( vec_string_type const & vs, string_view const where ) const
  {
    GC_assert(
      get_type() == GC_type::MAP,
      "{} bad data type, expect: {} but data stored is of type: {}",
      where,
      to_string( GC_type::MAP ),
      to_string( get_type() ) );
    string keys;
    for ( string_type const & s : vs )
    {
      keys += ' ';
      keys += s;
      if ( auto iv{ _m().find( s ) }; iv != _m().end() ) return s;
    }
    GC_assert( false, "{} cant find in keys: {}", where, keys );
    return "";
  }

  // -----------------------------------------------------------------------

  bool GenericContainer::get_if_exists( string_view const field, bool & value ) const
  {
    if ( get_type() != GC_type::MAP ) return false;
    auto const iv{ _m().find( string_type( field ) ) };
    if ( iv == _m().end() ) return false;
    if ( iv->second.get_type() != GC_type::BOOL ) return false;
    value = iv->second._b();
    return true;
  }

  bool GenericContainer::get_if_exists( vec_string_type const & fields, bool & value ) const
  {
    for ( string_view const field : fields )
    {
      if ( get_if_exists( field, value ) ) return true;
    }
    return false;
  }

  // -----------------------------------------------------------------------

  bool GenericContainer::get_if_exists( string_view const field, int_type & value ) const
  {
    if ( get_type() != GC_type::MAP ) return false;
    auto const iv{ _m().find( string_type( field ) ) };
    if ( iv == _m().end() ) return false;
    switch ( iv->second.get_type() )
    {
      case GC_type::BOOL: value = iv->second._b() ? 1 : 0; break;
      case GC_type::INTEGER: value = iv->second._i(); break;
      case GC_type::LONG: value = static_cast<int_type>( iv->second._l() ); break;
      case GC_type::REAL:
        if ( !GC_details::real_fits_integral<int_type>( iv->second._r() ) ) return false;
        value = static_cast<int_type>( iv->second._r() );
        break;
      case GC_type::COMPLEX:
        if ( !( GC_details::real_fits_integral<int_type>( iv->second._c().real() ) &&
                isZero0( iv->second._c().imag() ) ) )
          return false;
        value = static_cast<int_type>( iv->second._c().real() );
        break;
      case GC_type::NOTYPE:
      case GC_type::POINTER:
      case GC_type::STRING:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_BOOL:
      case GC_type::VEC_INTEGER:
      case GC_type::VEC_LONG:
      case GC_type::VEC_REAL:
      case GC_type::VEC_COMPLEX:
      case GC_type::VEC_STRING:
      case GC_type::MAT_INTEGER:
      case GC_type::MAT_LONG:
      case GC_type::MAT_REAL:
      case GC_type::MAT_COMPLEX:
      case GC_type::VECTOR:
      case GC_type::MAP: return false;
    }
    return true;
  }

  // -----------------------------------------------------------------------

  bool GenericContainer::get_if_exists( string_view const field, uint_type & value ) const
  {
    if ( get_type() != GC_type::MAP ) return false;
    auto const iv{ _m().find( string_type( field ) ) };
    if ( iv == _m().end() ) return false;
    switch ( iv->second.get_type() )
    {
      case GC_type::BOOL: value = iv->second._b() ? 1 : 0; break;
      case GC_type::INTEGER:
        if ( !std::in_range<uint_type>( iv->second._i() ) ) return false;
        value = static_cast<uint_type>( iv->second._i() );
        break;
      case GC_type::LONG:
        if ( !std::in_range<uint_type>( iv->second._l() ) ) return false;
        value = static_cast<uint_type>( iv->second._l() );
        break;
      case GC_type::REAL:
        if ( !GC_details::real_fits_integral<uint_type>( iv->second._r() ) ) return false;
        value = static_cast<uint_type>( iv->second._r() );
        break;
      case GC_type::COMPLEX:
        if ( !( GC_details::real_fits_integral<uint_type>( iv->second._c().real() ) &&
                isZero0( iv->second._c().imag() ) ) )
          return false;
        value = static_cast<uint_type>( iv->second._c().real() );
        break;
      case GC_type::NOTYPE:
      case GC_type::POINTER:
      case GC_type::STRING:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_BOOL:
      case GC_type::VEC_INTEGER:
      case GC_type::VEC_LONG:
      case GC_type::VEC_REAL:
      case GC_type::VEC_COMPLEX:
      case GC_type::VEC_STRING:
      case GC_type::MAT_INTEGER:
      case GC_type::MAT_LONG:
      case GC_type::MAT_REAL:
      case GC_type::MAT_COMPLEX:
      case GC_type::VECTOR:
      case GC_type::MAP: return false;
    }
    return true;
  }

  // -----------------------------------------------------------------------

  bool GenericContainer::get_if_exists( string_view const field, long_type & value ) const
  {
    if ( get_type() != GC_type::MAP ) return false;
    auto const iv{ _m().find( string_type( field ) ) };
    if ( iv == _m().end() ) return false;
    switch ( iv->second.get_type() )
    {
      case GC_type::BOOL: value = iv->second._b() ? 1 : 0; break;
      case GC_type::INTEGER: value = static_cast<long_type>( iv->second._i() ); break;
      case GC_type::LONG: value = iv->second._l(); break;
      case GC_type::REAL:
        if ( !GC_details::real_fits_integral<long_type>( iv->second._r() ) ) return false;
        value = static_cast<long_type>( iv->second._r() );
        break;
      case GC_type::COMPLEX:
        if ( !( GC_details::real_fits_integral<long_type>( iv->second._c().real() ) &&
                isZero0( iv->second._c().imag() ) ) )
          return false;
        value = static_cast<long_type>( iv->second._c().real() );
        break;
      case GC_type::NOTYPE:
      case GC_type::POINTER:
      case GC_type::STRING:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_BOOL:
      case GC_type::VEC_INTEGER:
      case GC_type::VEC_LONG:
      case GC_type::VEC_REAL:
      case GC_type::VEC_COMPLEX:
      case GC_type::VEC_STRING:
      case GC_type::MAT_INTEGER:
      case GC_type::MAT_LONG:
      case GC_type::MAT_REAL:
      case GC_type::MAT_COMPLEX:
      case GC_type::VECTOR:
      case GC_type::MAP: return false;
    }
    return true;
  }

  // -----------------------------------------------------------------------

  bool GenericContainer::get_if_exists( string_view const field, ulong_type & value ) const
  {
    if ( get_type() != GC_type::MAP ) return false;
    auto const iv{ _m().find( string_type( field ) ) };
    if ( iv == _m().end() ) return false;
    switch ( iv->second.get_type() )
    {
      case GC_type::BOOL: value = iv->second._b() ? 1 : 0; break;
      case GC_type::INTEGER:
        if ( !std::in_range<ulong_type>( iv->second._i() ) ) return false;
        value = static_cast<ulong_type>( iv->second._i() );
        break;
      case GC_type::LONG:
        if ( !std::in_range<ulong_type>( iv->second._l() ) ) return false;
        value = static_cast<ulong_type>( iv->second._l() );
        break;
      case GC_type::REAL:
        if ( !GC_details::real_fits_integral<ulong_type>( iv->second._r() ) ) return false;
        value = static_cast<ulong_type>( iv->second._r() );
        break;
      case GC_type::COMPLEX:
        if ( !( GC_details::real_fits_integral<ulong_type>( iv->second._c().real() ) &&
                isZero0( iv->second._c().imag() ) ) )
          return false;
        value = static_cast<ulong_type>( iv->second._c().real() );
        break;
      case GC_type::NOTYPE:
      case GC_type::POINTER:
      case GC_type::STRING:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_BOOL:
      case GC_type::VEC_INTEGER:
      case GC_type::VEC_LONG:
      case GC_type::VEC_REAL:
      case GC_type::VEC_COMPLEX:
      case GC_type::VEC_STRING:
      case GC_type::MAT_INTEGER:
      case GC_type::MAT_LONG:
      case GC_type::MAT_REAL:
      case GC_type::MAT_COMPLEX:
      case GC_type::VECTOR:
      case GC_type::MAP: return false;
    }
    return true;
  }

  // -----------------------------------------------------------------------

  bool GenericContainer::get_if_exists( string_view const field, real_type & value ) const
  {
    if ( get_type() != GC_type::MAP ) return false;
    auto const iv{ _m().find( string_type( field ) ) };
    if ( iv == _m().end() ) return false;
    switch ( iv->second.get_type() )
    {
      case GC_type::BOOL: value = static_cast<real_type>( iv->second._b() ? 1 : 0 ); break;
      case GC_type::INTEGER: value = static_cast<real_type>( iv->second._i() ); break;
      case GC_type::LONG: value = static_cast<real_type>( iv->second._l() ); break;
      case GC_type::REAL: value = iv->second._r(); break;
      case GC_type::COMPLEX:
        if ( !isZero0( iv->second._c().imag() ) ) return false;
        value = iv->second._c().real();
        break;
      case GC_type::NOTYPE:
      case GC_type::POINTER:
      case GC_type::STRING:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_BOOL:
      case GC_type::VEC_INTEGER:
      case GC_type::VEC_LONG:
      case GC_type::VEC_REAL:
      case GC_type::VEC_COMPLEX:
      case GC_type::VEC_STRING:
      case GC_type::MAT_INTEGER:
      case GC_type::MAT_LONG:
      case GC_type::MAT_REAL:
      case GC_type::MAT_COMPLEX:
      case GC_type::VECTOR:
      case GC_type::MAP: return false;
    }
    return true;
  }

  // -----------------------------------------------------------------------

  bool GenericContainer::get_if_exists( string_view const field, complex_type & value ) const
  {
    if ( get_type() != GC_type::MAP ) return false;
    auto const iv{ _m().find( string_type( field ) ) };
    if ( iv == _m().end() ) return false;
    switch ( iv->second.get_type() )
    {
      case GC_type::BOOL: value = complex_type( iv->second._b() ? 1 : 0, 0 ); break;
      case GC_type::INTEGER: value = complex_type( static_cast<real_type>( iv->second._i() ), 0 ); break;
      case GC_type::LONG: value = complex_type( static_cast<real_type>( iv->second._l() ), 0 ); break;
      case GC_type::REAL: value = complex_type( iv->second._r(), 0 ); break;
      case GC_type::COMPLEX: value = iv->second._c(); break;
      case GC_type::NOTYPE:
      case GC_type::POINTER:
      case GC_type::STRING:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_BOOL:
      case GC_type::VEC_INTEGER:
      case GC_type::VEC_LONG:
      case GC_type::VEC_REAL:
      case GC_type::VEC_COMPLEX:
      case GC_type::VEC_STRING:
      case GC_type::MAT_INTEGER:
      case GC_type::MAT_LONG:
      case GC_type::MAT_REAL:
      case GC_type::MAT_COMPLEX:
      case GC_type::VECTOR:
      case GC_type::MAP: return false;
    }
    return true;
  }

  // -----------------------------------------------------------------------

  bool GenericContainer::get_if_exists( string_view const field, string_type & value ) const
  {
    if ( get_type() != GC_type::MAP ) return false;
    auto const iv{ _m().find( string_type( field ) ) };
    if ( iv == _m().end() ) return false;
    if ( iv->second.get_type() != GC_type::STRING ) return false;
    value = iv->second._s();
    return true;
  }

  // --------------------------------------------------------------
  bool_type GenericContainer::get_bool_at( std::size_t const i )
  {
    if ( get_type() == GC_type::NOTYPE ) set_vec_bool();
    if ( get_type() == GC_type::VEC_BOOL )
    {
      CHECK_RESIZE( _v_b(), i );  // correct type, check size
      return _v_b()[i];
    }
    if ( get_type() != GC_type::VECTOR ) promote_to_vector();
    CHECK_RESIZE( _v(), i );
    return _v()[i].set_bool( false );
  }

  bool_type GenericContainer::get_bool_at( std::size_t const i, string_view const where ) const
  {
    ck( where, GC_type::VEC_BOOL );
    GC_assert( i < _v_b().size(), "{} get_bool_at( {} ) const, out of range", where, i );
    return _v_b()[i];
  }

  int_type & GenericContainer::get_int_at( std::size_t const i )
  {
    if ( get_type() == GC_type::NOTYPE )
      set_vec_int();
    else if ( get_type() == GC_type::BOOL || get_type() == GC_type::INTEGER || get_type() == GC_type::VEC_BOOL )
      promote_to_vec_int();
    if ( get_type() == GC_type::VEC_INTEGER )
    {
      CHECK_RESIZE( _v_i(), i );  // correct type, check size
      return _v_i()[i];
    }
    if ( get_type() != GC_type::VECTOR ) promote_to_vector();
    CHECK_RESIZE( _v(), i );
    if ( _v().size() <= i )
    {
      _v().resize( i + 1 );
      _v().back().set_int( 0 );
    }
    return _v()[i].get_int();
  }

  int_type const & GenericContainer::get_int_at( std::size_t const i, string_view const where ) const
  {
    ck( where, GC_type::VEC_INTEGER );
    GC_assert( i < _v_i().size(), "{} get_int_at( {} ) const, out of range", where, i );
    return _v_i()[i];
  }

  int_type & GenericContainer::get_int_at( std::size_t const i, std::size_t const j )
  {
    if ( get_type() == GC_type::NOTYPE )
      set_mat_int( i + 1, j + 1 );
    else if ( get_type() == GC_type::VEC_BOOL || get_type() == GC_type::VEC_INTEGER )
      promote_to_mat_int();
    GC_assert(
      GC_type::MAT_INTEGER == get_type(),
      "get_int_at( {}, {} ) bad data type\n"
      "expect: {}\n"
      "but data stored is of type: {}",
      i,
      j,
      to_string( GC_type::MAT_INTEGER ),
      to_string( get_type() ) );
    return _m_i()( i, j );
  }

  int_type const & GenericContainer::get_int_at( std::size_t const i, std::size_t const j, string_view const where )
    const
  {
    ck( where, GC_type::MAT_INTEGER );
    GC_assert(
      i < _m_i().num_rows() && j < _m_i().num_cols(),
      "{} get_int_at( {}, {} ) const, out of range",
      where,
      i,
      j );
    return _m_i()( i, j );
  }

  long_type & GenericContainer::get_long_at( std::size_t const i )
  {
    if ( get_type() == GC_type::NOTYPE )
      set_vec_long();
    else if (
      get_type() == GC_type::BOOL || get_type() == GC_type::INTEGER || get_type() == GC_type::LONG ||
      get_type() == GC_type::VEC_BOOL || get_type() == GC_type::VEC_INTEGER )
      promote_to_vec_long();
    if ( get_type() == GC_type::VEC_LONG )
    {
      CHECK_RESIZE( _v_l(), i );  // correct type, check size
      return _v_l()[i];
    }
    if ( get_type() != GC_type::VECTOR ) promote_to_vector();
    CHECK_RESIZE( _v(), i );
    return _v()[i].set_long( 0 );
  }

  long_type const & GenericContainer::get_long_at( std::size_t const i, string_view const where ) const
  {
    ck( where, GC_type::VEC_LONG );
    GC_assert( i < _v_l().size(), "{} get_long_at( {} ) const, out of range", where, i );
    return _v_l()[i];
  }

  long_type & GenericContainer::get_long_at( std::size_t const i, std::size_t const j )
  {
    if ( get_type() == GC_type::NOTYPE )
      set_mat_long( i + 1, j + 1 );
    else if ( get_type() == GC_type::VEC_BOOL || get_type() == GC_type::VEC_INTEGER || get_type() == GC_type::VEC_LONG )
      promote_to_mat_long();
    GC_assert(
      GC_type::MAT_LONG == get_type(),
      "get_long_at( {}, {} ) bad data type\n"
      "expect: {}\n"
      "but data stored is of type: {}",
      i,
      j,
      to_string( GC_type::MAT_LONG ),
      to_string( get_type() ) );
    return _m_l()( i, j );
  }

  long_type const & GenericContainer::get_long_at( std::size_t const i, std::size_t const j, string_view const where )
    const
  {
    ck( where, GC_type::MAT_LONG );
    GC_assert(
      i < _m_l().num_rows() && j < _m_l().num_cols(),
      "{} get_long_at( {}, {} ) const, out of range",
      where,
      i,
      j );
    return _m_l()( i, j );
  }

  real_type & GenericContainer::get_real_at( std::size_t const i )
  {
    if ( get_type() == GC_type::NOTYPE )
      set_vec_real();
    else if ( get_type() == GC_type::VEC_BOOL || get_type() == GC_type::VEC_INTEGER || get_type() == GC_type::VEC_LONG )
      promote_to_vec_real();
    if ( get_type() == GC_type::VEC_REAL )
    {
      CHECK_RESIZE( _v_r(), i );  // correct type, check size
      return _v_r()[i];
    }
    if ( get_type() != GC_type::VECTOR ) promote_to_vector();
    CHECK_RESIZE( _v(), i );
    return _v()[i].set_real( 0 );
  }

  real_type const & GenericContainer::get_real_at( std::size_t const i, string_view const where ) const
  {
    ck( where, GC_type::VEC_REAL );
    GC_assert( i < _v_r().size(), "{} get_real_at( {} ) const, out of range", where, i );
    return _v_r()[i];
  }

  real_type & GenericContainer::get_real_at( std::size_t const i, std::size_t const j )
  {
    if ( get_type() == GC_type::NOTYPE )
      set_mat_real( i + 1, j + 1 );
    else if (
      get_type() == GC_type::VEC_BOOL || get_type() == GC_type::VEC_INTEGER || get_type() == GC_type::VEC_LONG ||
      get_type() == GC_type::VEC_REAL )
      promote_to_mat_real();
    GC_assert(
      GC_type::MAT_REAL == get_type(),
      "get_real_at( {}, {} ) bad data type\n"
      "expect: {}\n"
      "but data stored is of type: {}",
      i,
      j,
      to_string( GC_type::MAT_REAL ),
      to_string( get_type() ) );
    return _m_r()( i, j );
  }

  real_type const & GenericContainer::get_real_at( std::size_t const i, std::size_t const j, string_view const where )
    const
  {
    ck( where, GC_type::MAT_REAL );
    GC_assert(
      i < _m_r().num_rows() && j < _m_r().num_cols(),
      "{} get_real_at( {}, {} ) const, out of range",
      where,
      i,
      j );
    return _m_r()( i, j );
  }

  complex_type & GenericContainer::get_complex_at( std::size_t const i )
  {
    if ( get_type() == GC_type::NOTYPE )
      set_vec_complex();
    else if (
      get_type() == GC_type::VEC_BOOL || get_type() == GC_type::VEC_INTEGER || get_type() == GC_type::VEC_LONG ||
      get_type() == GC_type::VEC_REAL )
      promote_to_vec_complex();
    if ( get_type() == GC_type::VEC_COMPLEX )
    {
      CHECK_RESIZE( _v_c(), i );  // correct type, check size
      return _v_c()[i];
    }
    if ( get_type() != GC_type::VECTOR ) promote_to_vector();
    CHECK_RESIZE( _v(), i );
    return _v()[i].set_complex( 0, 0 );
  }

  complex_type const & GenericContainer::get_complex_at( std::size_t const i, string_view const where ) const
  {
    ck( where, GC_type::VEC_COMPLEX );
    GC_assert( i < _v_c().size(), "{} get_complex_at( {} ) const, out of range", where, i );
    return _v_c()[i];
  }

  complex_type & GenericContainer::get_complex_at( std::size_t const i, std::size_t const j )
  {
    if ( get_type() == GC_type::NOTYPE )
      set_mat_complex( i + 1, j + 1 );
    else if (
      get_type() == GC_type::VEC_BOOL || get_type() == GC_type::VEC_INTEGER || get_type() == GC_type::VEC_LONG ||
      get_type() == GC_type::VEC_REAL || get_type() == GC_type::VEC_COMPLEX || get_type() == GC_type::MAT_REAL )
      promote_to_mat_complex();
    GC_assert(
      GC_type::MAT_COMPLEX == get_type(),
      "get_complex_at( {}, {} ) bad data type\n"
      "expect: {}\n"
      "but data stored is of type: {}",
      i,
      j,
      to_string( GC_type::MAT_COMPLEX ),
      to_string( get_type() ) );
    return _m_c()( i, j );
  }

  complex_type const & GenericContainer::get_complex_at(
    std::size_t const i,
    std::size_t const j,
    string_view const where ) const
  {
    ck( where, GC_type::MAT_COMPLEX );
    GC_assert(
      i < _m_c().num_rows() && j < _m_c().num_cols(),
      "{} get_complex_at( {}, {} ) const, out of range",
      where,
      i,
      j );
    return _m_c()( i, j );
  }

  string_type & GenericContainer::get_string_at( std::size_t const i )
  {
    if ( get_type() == GC_type::NOTYPE ) set_vec_string();
    if ( get_type() == GC_type::VEC_STRING )
    {
      CHECK_RESIZE( _v_s(), i );
      return _v_s()[i];
    }
    promote_to_vector();
    return ( *this )[i].set_string( "" );
  }

  string const & GenericContainer::get_string_at( std::size_t const i, string_view const where ) const
  {
    ck( where, GC_type::VEC_STRING );
    GC_assert( i < _v_s().size(), "{} get_string_at( {} ) const, out of range", where, i );
    return _v_s()[i];
  }

  GenericContainer & GenericContainer::get_gc_at( std::size_t const i )
  { return ( *this )[i]; }

  GenericContainer const & GenericContainer::get_gc_at( std::size_t const i, string_view const where ) const
  { return ( *this )( i, where ); }

  /*
  //   _        __
  //  (_)_ __  / _| ___
  //  | | '_ \| |_ / _ \
  //  | | | | |  _| (_) |
  //  |_|_| |_|_|  \___/
  */

  GenericContainer const & GenericContainer::info( ostream_type & stream ) const
  {
    switch ( get_type() )
    {
      case GC_type::NOTYPE: stream << "GenericContainer: No data stored\n"; break;
      case GC_type::POINTER: stream << "Generic pointer: " << _p() << '\n'; break;
      case GC_type::BOOL: stream << "Boolean: " << ( _b() ? "true" : "false" ) << '\n'; break;
      case GC_type::INTEGER: stream << "Integer: " << _i() << '\n'; break;
      case GC_type::LONG: stream << "Long: " << _l() << '\n'; break;
      case GC_type::REAL: stream << "Floating Point: " << _r() << '\n'; break;
      case GC_type::COMPLEX: stream << "Complex Floating Point: " << to_string( _c() ) << '\n'; break;
      case GC_type::STRING: stream << "String: " << _s().c_str() << '\n'; break;
      case GC_type::VEC_POINTER: stream << "Vector of generic pointer of size " << _v_p().size() << '\n'; break;
      case GC_type::VEC_BOOL: stream << "Vector of boolean of size " << _v_b().size() << '\n'; break;
      case GC_type::VEC_INTEGER: stream << "Vector of integer of size " << _v_i().size() << '\n'; break;
      case GC_type::VEC_LONG: stream << "Vector of long integer of size " << _v_l().size() << '\n'; break;
      case GC_type::VEC_REAL: stream << "Vector of floating point number of size " << _v_r().size() << '\n'; break;
      case GC_type::VEC_COMPLEX:
        stream << "Vector of complex floating point number of size " << _v_c().size() << '\n';
        break;
      case GC_type::VEC_STRING: stream << "Vector of string of size " << _v_s().size() << '\n'; break;
      case GC_type::MAT_INTEGER: _m_i().info( stream ); break;
      case GC_type::MAT_LONG: _m_l().info( stream ); break;
      case GC_type::MAT_REAL: _m_r().info( stream ); break;
      case GC_type::MAT_COMPLEX: _m_c().info( stream ); break;
      case GC_type::VECTOR: stream << "Vector of generic data type of size " << _v().size() << '\n'; break;
      case GC_type::MAP:
        stream << "Map\n";
        break;
        // default:
        //   stream << "Type N. " << _data_type << " not recognized\n";
        //   break;
    }
    return *this;
  }

  /*
  //                              _               __ __
  //    ___  _ __   ___ _ __ __ _| |_ ___  _ __  | _|_ |
  //   / _ \| '_ \ / _ \ '__/ _` | __/ _ \| '__| | | | |
  //  | (_) | |_) |  __/ | | (_| | || (_) | |    | | | |
  //   \___/| .__/ \___|_|  \__,_|\__\___/|_|    | | | |
  //        |_|                                  |__|__|
  */

  GenericContainer & GenericContainer::operator[]( std::size_t const i )
  {
    switch ( ck( GC_type::VECTOR ) )
    {
      case 0: break;          // data present
      default: set_vector();  // data must be allocated;
    }
    CHECK_RESIZE( _v(), i );
    return _v()[i];
  }

  GenericContainer const & GenericContainer::operator[]( std::size_t const i ) const
  {
    GC_assert(
      GC_type::VECTOR == get_type(),
      "operator [] integer argument = {}\n"
      "expect: {}\n"
      "but data stored is of type: {}",
      i,
      to_string( GC_type::VECTOR ),
      to_string( get_type() ) );
    GC_assert( i < _v().size(), "operator [] const, index {} out of range", i );
    return _v()[i];
  }

  /*
  //                              _                ____
  //    ___  _ __   ___ _ __ __ _| |_ ___  _ __   / /\ \
  //   / _ \| '_ \ / _ \ '__/ _` | __/ _ \| '__| | |  | |
  //  | (_) | |_) |  __/ | | (_| | || (_) | |    | |  | |
  //   \___/| .__/ \___|_|  \__,_|\__\___/|_|    | |  | |
  //        |_|                                   \_\/_/
  //
  */

  GenericContainer & GenericContainer::operator()( std::size_t const i, string_view const where )
  {
    ck( where, GC_type::VECTOR );
    GC_assert( i < _v().size(), "{} operator () const, index {} out of range", where, i );
    return _v()[i];
  }

  GenericContainer const & GenericContainer::operator()( std::size_t const i, string_view const where ) const
  {
    ck( where, GC_type::VECTOR );
    GC_assert( i < _v().size(), "{} operator () const, index {} out of range", where, i );
    return _v()[i];
  }

  GenericContainer & GenericContainer::operator()( string_view s, string_view const where )
  {
    string_type const key{ s };
    GC_assert(
      GC_type::MAP == get_type(),
      "{} operator (), with string argument ``{}''\n"
      "expect: {}\n"
      "but data stored is of type: {}",
      where,
      s,
      to_string( GC_type::MAP ),
      to_string( get_type() ) );
    auto iv{ _m().find( key ) };
    GC_assert( iv != _m().end(), "{} operator(): Cannot find key '{}'!\npossibile keys: {}", where, s, get_keys() );
    return iv->second;
  }

  GenericContainer const & GenericContainer::operator()( string_view s, string_view const where ) const
  {
    string_type const key{ s };
    GC_assert(
      GC_type::MAP == get_type(),
      "{}\n"
      "operator() const, with string argument ``{}''\n"
      "expect: {}\n"
      "but data stored is of type: {}",
      where,
      s,
      to_string( GC_type::MAP ),
      to_string( get_type() ) );
    auto iv{ _m().find( key ) };
    GC_assert(
      iv != _m().end(),
      "{} operator() const: Cannot find key '{}'!\npossibile keys: {}",
      where,
      s,
      get_keys() );
    return iv->second;
  }

  // ------------

  GenericContainer & GenericContainer::operator()( vec_string_type const & vs, string_view const where )
  {
    GC_assert(
      GC_type::MAP == get_type(),
      "{} operator (), with vector of string argument\n"
      "expect: {}\n"
      "but data stored is of type: {}",
      where,
      to_string( GC_type::MAP ),
      to_string( get_type() ) );
    map_type & m = _m();
    for ( string_type const & s : vs )
    {
      if ( auto iv{ m.find( s ) }; iv != m.end() ) return iv->second;
    }
    GC_assert( false, "{} operator(): Cannot find the key!", where );
    return *this;
  }

  GenericContainer const & GenericContainer::operator()( vec_string_type const & vs, string_view const where ) const
  {
    GC_assert(
      GC_type::MAP == get_type(),
      "{} operator (), with vector of string argument\n"
      "expect: {}\n"
      " but data stored is of type: {}",
      where,
      to_string( GC_type::MAP ),
      to_string( get_type() ) );
    map_type const & m = _m();
    for ( string_type const & s : vs )
    {
      if ( auto iv{ m.find( s ) }; iv != m.end() ) return iv->second;
    }
    GC_assert( false, "{}\noperator(): Cannot find the key!", where );
    return *this;
  }

  /*
  //                              _               __ __
  //    ___  _ __   ___ _ __ __ _| |_ ___  _ __  | _|_ |
  //   / _ \| '_ \ / _ \ '__/ _` | __/ _ \| '__| | | | |
  //  | (_) | |_) |  __/ | | (_| | || (_) | |    | | | |
  //   \___/| .__/ \___|_|  \__,_|\__\___/|_|    | | | |
  //        |_|                                  |__|__|
  */

  GenericContainer & GenericContainer::operator[]( string_view const s )
  {
    if ( ck( GC_type::MAP ) != 0 ) set_map();  // if not data present allocate!
    return _m()[string_type( s )];
  }

  GenericContainer const & GenericContainer::operator[]( string_view const s ) const
  {
    string_type const key{ s };
    GC_assert(
      GC_type::MAP == get_type(),
      "operator [] string argument ``{}''\n"
      "expect: {}\n"
      "but data stored is of type: {}",
      s,
      to_string( GC_type::MAP ),
      to_string( get_type() ) );
    auto const iv{ _m().find( key ) };
    if ( iv == _m().end() )
    {
      GC_assert( false, "operator [] const: Cannot find key '{}'!\npossible keys: {}", s, get_keys() );
    }
    return iv->second;
  }

  /*
  //              _       _
  //   _ __  _ __(_)_ __ | |_
  //  | '_ \| '__| | '_ \| __|
  //  | |_) | |  | | | | | |_
  //  | .__/|_|  |_|_| |_|\__|
  //  |_|
  */

  void GenericContainer::dump( ostream_type & stream, string_view const prefix, string_view const indent ) const
  {
    switch ( get_type() )
    {
      case GC_type::NOTYPE: stream << prefix << "null\n"; break;
      case GC_type::POINTER:
        stream << prefix << std::hex << std::showbase << reinterpret_cast<uintptr_t>( _p() ) << '\n';
        break;
      case GC_type::BOOL: stream << prefix << ( _b() ? "true" : "false" ) << '\n'; break;
      case GC_type::INTEGER: stream << prefix << _i() << '\n'; break;
      case GC_type::LONG: stream << prefix << _l() << '\n'; break;
      case GC_type::REAL: stream << prefix << _r() << '\n'; break;
      case GC_type::COMPLEX: stream << prefix << to_string( _c() ) << '\n'; break;
      case GC_type::STRING: stream << prefix << "\"" << _s() << "\"\n"; break;
      case GC_type::VEC_POINTER:
      {
        vec_pointer_type const & v{ _v_p() };
        for ( vec_pointer_type::size_type i{ 0 }; i < v.size(); ++i )
          stream << prefix << "vec_pointer(" << i << "): " << std::hex << std::showbase
                 << reinterpret_cast<uintptr_t>( v[i] ) << '\n';
      }
      break;
      case GC_type::VEC_BOOL:
      {
        vec_bool_type const & v{ _v_b() };
        stream << prefix << v << '\n';
      }
      break;
      case GC_type::VEC_INTEGER:
      {
        vec_int_type const & v{ _v_i() };
        stream << prefix << v << '\n';
      }
      break;
      case GC_type::VEC_LONG:
      {
        vec_long_type const & v{ _v_l() };
        stream << prefix << v << '\n';
      }
      break;
      case GC_type::VEC_REAL:
      {
        vec_real_type const & v{ _v_r() };
        stream << prefix << v << '\n';
      }
      break;
      case GC_type::VEC_COMPLEX:
      {
        vec_complex_type const & v{ _v_c() };
        stream << prefix << v << '\n';
      }
      break;
      case GC_type::MAT_INTEGER:
      {
        mat_int_type const & m{ _m_i() };
        stream << m;
      }
      break;
      case GC_type::MAT_LONG:
      {
        mat_long_type const & m{ _m_l() };
        stream << m;
      }
      break;
      case GC_type::MAT_REAL:
      {
        mat_real_type const & m{ _m_r() };
        stream << m;
      }
      break;
      case GC_type::MAT_COMPLEX:
      {
        mat_complex_type const & m{ _m_c() };
        stream << m;
      }
      break;
      case GC_type::VEC_STRING:
      {
        vec_string_type const & v{ _v_s() };
        for ( vec_string_type::size_type i{ 0 }; i < v.size(); ++i ) stream << prefix << i << ": \"" << v[i] << "\"\n";
      }
      break;

      case GC_type::VECTOR:
      {
        vector_type const & v{ _v() };
        string              prefix2{ prefix };
        prefix2 += indent;
        for ( vector_type::size_type i{ 0 }; i < v.size(); ++i )
        {
          GenericContainer const & vi{ v[i] };
          if ( vi.simple_data() || ( vi.simple_vec_data() && vi.get_num_elements() <= 10 ) )
          {
            stream << prefix << i << ": ";
            vi.dump( stream, "" );
          }
          else
          {
            stream << prefix << i << ":\n";
            vi.dump( stream, prefix2 );
          }
        }
      }
      break;
      case GC_type::MAP:
      {
        map_type const & m{ _m() };
        string           prefix2{ prefix };
        prefix2 += indent;
        for ( const auto & [fst, snd] : m )
        {
          // check formatting using pcre
          // num+"@"+"underline character"
          // Try to find the regex in aLineToMatch, and report results.
          string_type matches[4];
          if ( int const pcreExecRet{ pcre_for_GC.exec( fst, matches ) }; pcreExecRet == 4 )
          {
            string_type header = matches[3];  // header
            // found formatting
            if ( snd.simple_data() )
            {
              stream << prefix << header << ": ";
              snd.dump( stream, "" );
            }
            else
            {
              if ( matches[1].length() > 1 ) stream << '\n';  // double ## --> add nel line
              stream << prefix << header;
              if ( matches[2].length() > 0 )
              {
                stream << '\n' << prefix;
                char const  fmt{ matches[2][0] };  // underline char
                std::size_t m3 = header.length();
                while ( m3-- > 0 ) stream << fmt;  // underline header
              }
              else
              {
                stream << ':';
              }
              stream << '\n';
              snd.dump( stream, prefix2 );
            }
          }
          else
          {
            string_type header = pcreExecRet == 3 ? matches[3] : fst;
            if ( snd.simple_data() )
            {
              stream << prefix << header << ": ";
              snd.dump( stream, "" );
            }
            else
            {
              stream << prefix << header << ":\n";
              snd.dump( stream, prefix2 );
            }
          }
        }
      }
      break;

        // default:
        //   GC_assert( false, "Error, print(...) unknown type!\n");
        //   break;
    }
  }


  void GenericContainer::print_content_types(
    ostream_type &    stream,
    string_view const prefix,
    string_view const indent ) const
  {
    switch ( get_type() )
    {
      case GC_type::NOTYPE: stream << prefix << "Empty!\n"; break;
      case GC_type::POINTER: stream << prefix << "(*void)\n"; break;
      case GC_type::BOOL: stream << prefix << "bool\n"; break;
      case GC_type::INTEGER: stream << prefix << "int\n"; break;
      case GC_type::LONG: stream << prefix << "long int\n"; break;
      case GC_type::REAL: stream << prefix << "double\n"; break;
      case GC_type::COMPLEX: stream << prefix << "complex\n"; break;
      case GC_type::STRING: stream << prefix << "string\n"; break;
      case GC_type::VEC_POINTER:
      {
        vec_pointer_type const & v{ _v_p() };
        stream << "vector of pointer[" << v.size() << "]\n";
      }
      break;
      case GC_type::VEC_BOOL:
      {
        vec_bool_type const & v{ _v_b() };
        stream << "vector of bool[" << v.size() << "]\n";
      }
      break;
      case GC_type::VEC_INTEGER:
      {
        vec_int_type const & v{ _v_i() };
        stream << "vector of int[" << v.size() << "]\n";
      }
      break;
      case GC_type::VEC_LONG:
      {
        vec_long_type const & v{ _v_l() };
        stream << "vector of long[" << v.size() << "]\n";
      }
      break;
      case GC_type::VEC_REAL:
      {
        vec_real_type const & v{ _v_r() };
        stream << "vector of double[" << v.size() << "]\n";
      }
      break;
      case GC_type::VEC_COMPLEX:
      {
        vec_complex_type const & v{ _v_c() };
        stream << "vector of complex[" << v.size() << "]\n";
      }
      break;
      case GC_type::MAT_INTEGER:
      {
        mat_int_type const & m{ _m_i() };
        stream << "matrix of int[" << m.num_rows() << "," << m.num_cols() << "]\n";
      }
      break;
      case GC_type::MAT_LONG:
      {
        mat_long_type const & m{ _m_l() };
        stream << "matrix of long[" << m.num_rows() << "," << m.num_cols() << "]\n";
      }
      break;
      case GC_type::MAT_REAL:
      {
        mat_real_type const & m{ _m_r() };
        stream << "matrix of double[" << m.num_rows() << "," << m.num_cols() << "]\n";
      }
      break;
      case GC_type::MAT_COMPLEX:
      {
        mat_complex_type const & m{ _m_c() };
        stream << "matrix of complex[" << m.num_rows() << "," << m.num_cols() << "]\n";
      }
      break;
      case GC_type::VEC_STRING:
      {
        vec_string_type const & v{ _v_s() };
        stream << "vector of string[" << v.size() << "]\n";
      }
      break;
      case GC_type::VECTOR:
      {
        vector_type const & v{ _v() };
        string              prefix2{ prefix };
        prefix2 += indent;
        for ( vector_type::size_type i{ 0 }; i < v.size(); ++i )
        {
          GenericContainer const & vi{ v[i] };
          if ( vi.simple_data() || vi.simple_vec_data() )
          {
            stream << prefix << i << ": ";
            vi.print_content_types( stream, "" );
          }
          else
          {
            stream << prefix << i << ":\n";
            vi.print_content_types( stream, prefix2, indent );
          }
        }
      }
      break;
      case GC_type::MAP:
      {
        map_type const & m{ _m() };
        string           prefix2{ prefix };
        prefix2 += indent;
        for ( const auto & [fst, snd] : m )
        {
          // check formatting using pcre
          // num+"@"+"underline character"
          // Try to find the regex in aLineToMatch, and report results.
          string_type matches[4];
          if ( int const pcreExecRet{ pcre_for_GC.exec( fst, matches ) }; pcreExecRet == 4 )
          {
            string_type header = matches[3];  // header
            // found formatting
            if ( snd.simple_data() || snd.simple_vec_data() )
            {
              stream << prefix << header << ": ";
              snd.print_content_types( stream, "" );
            }
            else
            {
              if ( matches[1].length() > 1 ) stream << '\n';  // double ## --> add nel line
              stream << prefix << header;
              if ( !matches[2].empty() )
              {
                stream << '\n' << prefix;
                char const  fmt{ matches[2][0] };  // underline char
                std::size_t m3 = header.length();
                while ( m3-- > 0 ) stream << fmt;  // underline header
              }
              else
              {
                stream << ':';
              }
              stream << '\n';
              snd.print_content_types( stream, prefix2, indent );
            }
          }
          else
          {
            string_type header = pcreExecRet == 3 ? matches[3] : fst;
            if ( snd.simple_data() || snd.simple_vec_data() )
            {
              stream << prefix << header << ": ";
              snd.print_content_types( stream, "" );
            }
            else
            {
              stream << prefix << header << ":\n";
              snd.print_content_types( stream, prefix2, indent );
            }
          }
        }
      }
      break;

        // default:
        //   GC_assert( false, "Error, print(...) unknown type!\n");
        //   break;
    }
  }

  /*
  //   _
  //  | |_ ___      __ _  ___
  //  | __/ _ \    / _` |/ __|
  //  | || (_) |  | (_| | (__
  //   \__\___/____\__, |\___|
  //         |_____|___/
  */

  void GenericContainer::to_gc( GenericContainer & gc ) const
  {
    switch ( get_type() )
    {
      case GC_type::NOTYPE: gc.clear(); break;
      case GC_type::BOOL: gc = _b(); break;
      case GC_type::INTEGER: gc = _i(); break;
      case GC_type::LONG: gc = _l(); break;
      case GC_type::REAL: gc = _r(); break;
      case GC_type::COMPLEX: gc = _c(); break;
      case GC_type::STRING: gc = _s(); break;
      case GC_type::VEC_BOOL: gc = _v_b(); break;
      case GC_type::VEC_INTEGER: gc = _v_i(); break;
      case GC_type::VEC_LONG: gc = _v_l(); break;
      case GC_type::VEC_REAL: gc = _v_r(); break;
      case GC_type::VEC_STRING: gc = _v_s(); break;

      case GC_type::VECTOR:
        gc.set_vector();
        {
          vector_type const & v{ _v() };
          vector_type &       vv = gc.set_vector( static_cast<std::size_t>( v.size() ) );
          for ( vector_type::size_type i{ 0 }; i < v.size(); ++i ) v[i].to_gc( vv[i] );
        }
        break;
      case GC_type::MAP:
        gc.set_map();
        {
          map_type const & m{ _m() };
          for ( const auto & [fst, snd] : m ) snd.to_gc( gc[fst] );
        }
        break;
      case GC_type::MAT_INTEGER: gc = _m_i(); break;
      case GC_type::MAT_LONG: gc = _m_l(); break;
      case GC_type::MAT_REAL: gc = _m_r(); break;
      case GC_type::VEC_COMPLEX: gc = _v_c(); break;
      case GC_type::MAT_COMPLEX: gc = _m_c(); break;
      case GC_type::POINTER: gc = this->get_pointer<void *>(); break;
      case GC_type::VEC_POINTER:
      {
        vec_pointer_type const & v  = this->get_vec_pointer();
        vec_pointer_type &       vv = gc.set_vec_pointer( v.size() );
        for ( vec_pointer_type::size_type i{ 0 }; i < v.size(); ++i ) vv[i] = v[i];
      }
      break;
    }
  }

  /*
  //    __
  //   / _|_ __ ___  _ __ ___       __ _  ___
  //  | |_| '__/ _ \| '_ ` _ \     / _` |/ __|
  //  |  _| | | (_) | | | | | |   | (_| | (__
  //  |_| |_|  \___/|_| |_| |_|____\__, |\___|
  //                         |_____|___/
  */

  void GenericContainer::from_gc( GenericContainer const & gc )
  {
    switch ( gc.get_type() )
    {
      case GC_type::NOTYPE: this->clear(); break;
      case GC_type::BOOL: this->set_bool( gc._b() ); break;
      case GC_type::INTEGER: this->set_int( gc._i() ); break;
      case GC_type::LONG: this->set_long( gc._l() ); break;
      case GC_type::REAL: this->set_real( gc._r() ); break;
      case GC_type::COMPLEX: this->set_complex( gc._c() ); break;
      case GC_type::STRING: this->set_string( gc._s() ); break;
      case GC_type::VEC_BOOL: this->set_vec_bool( gc._v_b() ); break;
      case GC_type::VEC_INTEGER: this->set_vec_int( gc._v_i() ); break;
      case GC_type::VEC_LONG: this->set_vec_long( gc._v_l() ); break;
      case GC_type::VEC_REAL: this->set_vec_real( gc._v_r() ); break;
      case GC_type::VEC_STRING: this->set_vec_string( gc._v_s() ); break;

      case GC_type::VECTOR:
        this->set_vector();
        {
          vector_type const & v{ gc._v() };
          vector_type &       vv = this->set_vector( static_cast<std::size_t>( v.size() ) );
          for ( vector_type::size_type i{ 0 }; i < v.size(); ++i ) vv[i].from_gc( v[i] );
        }
        break;
      case GC_type::MAP:
        this->set_map();
        {
          map_type const & m{ gc._m() };
          for ( const auto & [fst, snd] : m ) ( *this )[fst].from_gc( snd );
        }
        break;
      case GC_type::MAT_INTEGER: this->set_mat_int( gc._m_i() ); break;
      case GC_type::MAT_LONG: this->set_mat_long( gc._m_l() ); break;
      case GC_type::MAT_REAL: this->set_mat_real( gc._m_r() ); break;
      case GC_type::VEC_COMPLEX: this->set_vec_complex( gc._v_c() ); break;
      case GC_type::MAT_COMPLEX: this->set_mat_complex( gc._m_c() ); break;
      case GC_type::POINTER: this->set_pointer( gc.get_pointer<void *>() ); break;
      case GC_type::VEC_POINTER:
      {
        vec_pointer_type const & v  = gc.get_vec_pointer();
        vec_pointer_type &       vv = this->set_vec_pointer( v.size() );
        for ( vec_pointer_type::size_type i{ 0 }; i < v.size(); ++i ) vv[i] = v[i];
      }
      break;
    }
  }

  //
  // -------------------------------------------------------------------
  //

  void GenericContainer::merge( GenericContainer const & gc, string_view const where )
  {
    if ( gc.get_type() == GC_type::NOTYPE ) return;
    GC_assert(
      gc.get_type() == GC_type::MAP,
      "{} in merge data expected to be of type: {}\n"
      "but data stored is of type: {}",
      where,
      to_string( GC_type::MAP ),
      gc.get_type_name() );
    if ( get_type() == GC_type::NOTYPE ) this->set_map();
    ck( where, GC_type::MAP );
    {
      map_type const & m{ gc.get_map() };
      for ( const auto & [fst, snd] : m ) ( *this )[fst].from_gc( snd );
    }
  }
}  // namespace GC_namespace

//
// eof: GenericContainer.cc
//
