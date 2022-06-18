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
/// file: Numbers.hxx
///

namespace Utils {

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  using std::numeric_limits;
  using std::string;
  using std::sqrt;
  using std::fpclassify;
  using std::floor;
  #endif

  /*
  //    ____                _              _
  //   / ___|___  _ __  ___| |_ __ _ _ __ | |_ ___
  //  | |   / _ \| '_ \/ __| __/ _` | '_ \| __/ __|
  //  | |__| (_) | | | \__ \ || (_| | | | | |_\__ \
  //   \____\___/|_| |_|___/\__\__,_|_| |_|\__|___/
  */

  /// Not a number constant
  template <typename T> T NaN();
  template <typename T> T Inf();
  /// machine epsilon
  template <typename T> T machine_eps();

  /// square root of machine epsilon
  template <typename T> T sqrt_machine_eps();

  /// maximum representable value
  template <typename T> T maximum_value();

  /// minimum representable value
  template <typename T> T minimum_value();

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  template <> inline float
  NaN() { return numeric_limits<float>::quiet_NaN(); }

  template <> inline double
  NaN() { return numeric_limits<double>::quiet_NaN(); }

  template <> inline float
  Inf() { return numeric_limits<float>::infinity(); }

  template <> inline double
  Inf() { return numeric_limits<double>::infinity(); }

  template <> inline float
  machine_eps() { return numeric_limits<float>::epsilon(); }

  template <> inline double
  machine_eps() { return numeric_limits<double>::epsilon(); }

  template <> inline float
  sqrt_machine_eps() { return sqrt(numeric_limits<float>::epsilon()); }

  template <> inline double
  sqrt_machine_eps() { return sqrt(numeric_limits<double>::epsilon()); }

  template <> inline float
  maximum_value() { return sqrt(numeric_limits<float>::max()); }

  template <> inline double
  maximum_value() { return sqrt(numeric_limits<double>::max()); }

  template <> inline float
  minimum_value() { return sqrt(numeric_limits<float>::min()); }

  template <> inline double
  minimum_value() { return sqrt(numeric_limits<double>::min()); }
  #endif

  static inline bool is_zero( double x )     { return FP_ZERO == fpclassify(x); }
  static inline bool is_zero( float x )      { return FP_ZERO == fpclassify(x); }
  static inline bool is_NaN( double x )      { return std::isnan(x); }
  static inline bool is_NaN( float x )       { return std::isnan(x); }
  static inline bool is_infinite( double x ) { return std::isinf(x); }
  static inline bool is_infinite( float x )  { return std::isinf(x); }
  static inline bool is_normal( double x )   { return std::isnormal(x); }
  static inline bool is_normal( float x )    { return std::isnormal(x); }
  static inline bool is_finite( double x )   { return std::isfinite(x); }
  static inline bool is_finite( float x )    { return std::isfinite(x); }
  static inline bool is_integer( double x )  { return is_zero( x-floor(x) ); }
  static inline bool is_integer( float x )   { return is_zero( x-floor(x) ); }
  static inline bool is_unsigned( double x ) { return is_integer(x) && x >= 0; }
  static inline bool is_unsigned( float x )  { return is_integer(x) && x >= 0; }

  //============================================================================

  bool found_NaN( double const * pv, int DIM );
  bool found_NaN( float const * pv, int DIM );

  void
  check_NaN(
    double const * pv,
    char   const * v_name,
    int            DIM,
    int            line,
    char   const * file
  );

  void
  check_NaN(
    float const * pv,
    char  const * v_name,
    int           DIM,
    int           line,
    char  const * file
  );

  //============================================================================

  //!
  //! `m_e` the value of \f$ e \f$.
  //!
  static double const m_e = 2.718281828459045235360287471352662497757;

  //!
  //! `m_pi` the value of \f$ \pi \f$.
  //!
  static double const m_pi = 3.141592653589793238462643383279502884197;

  //!
  //! `m_2pi` the value of \f$ 2\pi \f$.
  //!
  static double const m_2pi = 6.283185307179586476925286766559005768394;

  //!
  //! `m_pi_2` the value of \f$ \pi/2 \f$.
  //!
  static double const m_pi_2 = 1.570796326794896619231321691639751442098;

  //!
  //! `m_pi_4` the value of \f$ \pi/4 \f$.
  //!
  static double const m_pi_4 = 0.7853981633974483096156608458198757210492;

  //!
  //! `m_1_pi` the value of \f$ 1/\pi \f$.
  //!
  static double const m_1_pi = 0.3183098861837906715377675267450287240689;

  //!
  //! `m_2_pi` the value of \f$ 2/\pi \f$.
  //!
  static double const m_2_pi = 0.6366197723675813430755350534900574481378;

  //!
  //! `m_sqrtpi` the value of \f$ \sqrt{\pi} \f$.
  //!
  static double const m_sqrtpi = 1.772453850905516027298167483341145182798;

  //!
  //! `m_2_sqrtpi` the value of \f$ 2/\sqrt{\pi} \f$.
  //!
  static double const m_2_sqrtpi = 1.128379167095512573896158903121545171688;

  //!
  //! `m_sqrt2` the value of \f$ \sqrt{2} \f$.
  //!
  static double const m_sqrt2 = 1.414213562373095048801688724209698078570;

  //!
  //! `m_1_sqrt2` the value of \f$ 1/\sqrt{2} \f$.
  //!
  static double const m_1_sqrt2 = 0.7071067811865475244008443621048490392850;

  #ifdef UTILS_OLD_CAMELCASE
  template <typename T> inline T machineEps() { return machine_eps<T>(); }
  template <typename T> inline T sqrtMachineEps() { return sqrt_machine_eps<T>(); }
  template <typename T> inline T maximumValue() { return maximum_value<T>(); }
  template <typename T> inline T minimumValue() { return minimum_value<T>(); }

  static inline bool isZero( double x )     { return is_zero(x); }
  static inline bool isZero( float x )      { return is_zero(x); }
  static inline bool isInfinite( double x ) { return is_infinite(x); }
  static inline bool isInfinite( float x )  { return is_infinite(x); }
  static inline bool isNaN( double x )      { return is_NaN(x); }
  static inline bool isNaN( float x )       { return is_NaN(x); }
  static inline bool isFinite( double x )   { return is_finite(x); }
  static inline bool isFinite( float x )    { return is_finite(x); }
  static inline bool isRegular( double x )  { return is_finite(x); }
  static inline bool isRegular( float x )   { return is_finite(x); }
  static inline bool isInteger( double x )  { return is_integer(x); }
  static inline bool isInteger( float x )   { return is_integer(x); }
  static inline bool isUnsigned( double x ) { return is_unsigned(x); }
  static inline bool isUnsigned( float x )  { return is_unsigned(x); }

  static inline bool
  foundNaN( double const * pv, int DIM )
  { return found_NaN( pv, DIM ); }

  static inline bool
  foundNaN( float const * pv, int DIM )
  { return found_NaN( pv, DIM ); }

  static inline void
  checkNaN(
    double const * pv,
    char   const * v_name,
    int            DIM,
    int            line,
    char   const * file
  ) {
    check_NaN( pv, v_name, DIM, line, file );
  }

  static inline void
  checkNaN(
    float const * pv,
    char  const * v_name,
    int           DIM,
    int           line,
    char  const * file
  ) {
    check_NaN( pv, v_name, DIM, line, file );
  }

  #endif
}

///
/// eof: Number.hxx
///
