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
// file: Numbers.hxx
//

namespace Utils
{

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  using std::floor;
  using std::fpclassify;
  using std::numeric_limits;
  using std::sqrt;
  using std::string;
#endif

  /*
  //    ____                _              _
  //   / ___|___  _ __  ___| |_ __ _ _ __ | |_ ___
  //  | |   / _ \| '_ \/ __| __/ _` | '_ \| __/ __|
  //  | |__| (_) | | | \__ \ || (_| | | | | |_\__ \
  //   \____\___/|_| |_|___/\__\__,_|_| |_|\__|___/
  */

  //! Returns a Not-a-Number (NaN) constant for the specified type.
  //!
  //! \tparam T The type of the NaN constant. Supported types include float and
  //! double.
  //! \return A NaN constant of type T.
  template <typename T>
  T NaN();

  //! Returns an infinity constant for the specified type.
  //!
  //! \tparam T The type of the infinity constant. Supported types include float
  //! and double.
  //! \return An infinity constant of type T.
  template <typename T>
  T Inf();

  //! Returns the machine epsilon for the specified type.
  //!
  //! \tparam T The type of the machine epsilon. Supported types include float
  //! and double.
  //! \return The machine epsilon constant of type T.
  template <typename T>
  T machine_eps();

  //! Returns the square root of the machine epsilon for the specified type.
  //!
  //! \tparam T The type for which to compute the square root of machine
  //! epsilon. Supported types include float and double.
  //! \return The square root of the machine epsilon of type T.
  template <typename T>
  T sqrt_machine_eps();

  //! Returns the maximum representable value for the specified type.
  //!
  //! \tparam T The type for which to compute the maximum value. Supported types
  //! include float and double.
  //! \return The maximum representable value of type T.
  template <typename T>
  T maximum_value();

  //! Returns the minimum representable value for the specified type.
  //!
  //! \tparam T The type for which to compute the minimum value. Supported types
  //! include float and double.
  //! \return The minimum representable value of type T.
  template <typename T>
  T minimum_value();

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  template <>
  inline float
  NaN()
  {
    return numeric_limits<float>::quiet_NaN();
  }
  template <>
  inline double
  NaN()
  {
    return numeric_limits<double>::quiet_NaN();
  }

  template <>
  inline float
  Inf()
  {
    return numeric_limits<float>::infinity();
  }
  template <>
  inline double
  Inf()
  {
    return numeric_limits<double>::infinity();
  }

  template <>
  inline float
  machine_eps()
  {
    return numeric_limits<float>::epsilon();
  }
  template <>
  inline double
  machine_eps()
  {
    return numeric_limits<double>::epsilon();
  }

  template <>
  inline float
  sqrt_machine_eps()
  {
    return sqrt( numeric_limits<float>::epsilon() );
  }
  template <>
  inline double
  sqrt_machine_eps()
  {
    return sqrt( numeric_limits<double>::epsilon() );
  }

  template <>
  inline float
  maximum_value()
  {
    return sqrt( numeric_limits<float>::max() );
  }
  template <>
  inline double
  maximum_value()
  {
    return sqrt( numeric_limits<double>::max() );
  }

  template <>
  inline float
  minimum_value()
  {
    return sqrt( numeric_limits<float>::min() );
  }
  template <>
  inline double
  minimum_value()
  {
    return sqrt( numeric_limits<double>::min() );
  }
#endif

  //! Checks if a double value is zero.
  //!
  //! \param x The double value to check.
  //! \return True if x is zero, otherwise false.
  static inline bool
  is_zero( double x )
  {
    return FP_ZERO == fpclassify( x );
  }

  //! Checks if a float value is zero.
  //!
  //! \param x The float value to check.
  //! \return True if x is zero, otherwise false.
  static inline bool
  is_zero( float x )
  {
    return FP_ZERO == fpclassify( x );
  }

  //! Checks if a double value is NaN (Not-a-Number).
  //!
  //! \param x The double value to check.
  //! \return True if x is NaN, otherwise false.
  static inline bool
  is_NaN( double x )
  {
    return std::isnan( x );
  }

  //! Checks if a float value is NaN (Not-a-Number).
  //!
  //! \param x The float value to check.
  //! \return True if x is NaN, otherwise false.
  static inline bool
  is_NaN( float x )
  {
    return std::isnan( x );
  }

  //! Checks if a double value is infinite.
  //!
  //! \param x The double value to check.
  //! \return True if x is infinite, otherwise false.
  static inline bool
  is_infinite( double x )
  {
    return std::isinf( x );
  }

  //! Checks if a float value is infinite.
  //!
  //! \param x The float value to check.
  //! \return True if x is infinite, otherwise false.
  static inline bool
  is_infinite( float x )
  {
    return std::isinf( x );
  }

  //! Checks if a double value is normal (not NaN or infinite).
  //!
  //! \param x The double value to check.
  //! \return True if x is normal, otherwise false.
  static inline bool
  is_normal( double x )
  {
    return std::isnormal( x );
  }

  //! Checks if a float value is normal (not NaN or infinite).
  //!
  //! \param x The float value to check.
  //! \return True if x is normal, otherwise false.
  static inline bool
  is_normal( float x )
  {
    return std::isnormal( x );
  }

  //! Checks if a double value is finite (not NaN or infinite).
  //!
  //! \param x The double value to check.
  //! \return True if x is finite, otherwise false.
  static inline bool
  is_finite( double x )
  {
    return std::isfinite( x );
  }

  //! Checks if a float value is finite (not NaN or infinite).
  //!
  //! \param x The float value to check.
  //! \return True if x is finite, otherwise false.
  static inline bool
  is_finite( float x )
  {
    return std::isfinite( x );
  }

  //! Checks if a double value is an integer.
  //!
  //! \param x The double value to check.
  //! \return True if x is an integer, otherwise false.
  static inline bool
  is_integer( double x )
  {
    return is_zero( x - floor( x ) );
  }

  //! Checks if a float value is an integer.
  //!
  //! \param x The float value to check.
  //! \return True if x is an integer, otherwise false.
  static inline bool
  is_integer( float x )
  {
    return is_zero( x - floor( x ) );
  }

  //! Checks if a double value is unsigned (non-negative integer).
  //!
  //! \param x The double value to check.
  //! \return True if x is unsigned, otherwise false.
  static inline bool
  is_unsigned( double x )
  {
    return is_integer( x ) && x >= 0;
  }

  //! Checks if a float value is unsigned (non-negative integer).
  //!
  //! \param x The float value to check.
  //! \return True if x is unsigned, otherwise false.
  static inline bool
  is_unsigned( float x )
  {
    return is_integer( x ) && x >= 0;
  }

  //============================================================================

  //! Checks if a NaN value is found in an array of doubles.
  //!
  //! \param pv Pointer to the array of double values.
  //! \param DIM The dimension (size) of the array.
  //! \return True if a NaN value is found, otherwise false.
  bool found_NaN( double const * pv, int DIM );

  //! Checks if a NaN value is found in an array of floats.
  //!
  //! \param pv Pointer to the array of float values.
  //! \param DIM The dimension (size) of the array.
  //! \return True if a NaN value is found, otherwise false.
  bool found_NaN( float const * pv, int DIM );

  //! Checks for NaN values in an array of doubles and logs an error if found.
  //!
  //! \param pv Pointer to the array of double values.
  //! \param v_name The name of the variable for logging.
  //! \param DIM The dimension (size) of the array.
  //! \param line The line number where the check is performed.
  //! \param file The name of the file where the check is performed.
  void check_NaN( double const * pv, string_view v_name, int DIM, int line, string_view file );

  //! Checks for NaN values in an array of floats and logs an error if found.
  //!
  //! \param pv Pointer to the array of float values.
  //! \param v_name The name of the variable for logging.
  //! \param DIM The dimension (size) of the array.
  //! \param line The line number where the check is performed.
  //! \param file The name of the file where the check is performed.
  void check_NaN( float const * pv, string_view v_name, int DIM, int line, string_view file );

  //============================================================================

  //! The value of \f$ e \f$ (Euler's number).
  static double const m_e = 2.718281828459045235360287471352662497757;

  //! The value of \f$ \pi \f$ (Pi).
  static double const m_pi = 3.141592653589793238462643383279502884197;

  //! The value of \f$ 2\pi \f$ (Two Pi).
  static double const m_2pi = 6.283185307179586476925286766559005768394;

  //! The value of \f$ \pi/2 \f$ (Pi divided by two).
  static double const m_pi_2 = 1.570796326794896619231321691639751442098;

  //! The value of \f$ \pi/4 \f$ (Pi divided by four).
  static double const m_pi_4 = 0.7853981633974483096156608458198757210492;

  //! The value of \f$ 1/\pi \f$ (One divided by Pi).
  static double const m_1_pi = 0.3183098861837906715377675267450287240689;

  //! The value of \f$ 2/\pi \f$ (Two divided by Pi).
  static double const m_2_pi = 0.6366197723675813430755350534900574481378;

  //! The value of \f$ \sqrt{\pi} \f$ (Square root of Pi).
  static double const m_sqrtpi = 1.772453850905516027298167483341145182798;

  //! The value of \f$ 2/\sqrt{\pi} \f$ (Two divided by the square root of Pi).
  static double const m_2_sqrtpi = 1.128379167095512573896158903121545171688;

  //! The value of \f$ \sqrt{2} \f$ (Square root of Two).
  static double const m_sqrt2 = 1.414213562373095048801688724209698078570;

  //! The value of \f$ 1/\sqrt{2} \f$ (One divided by the square root of Two).
  static double const m_1_sqrt2 = 0.7071067811865475244008443621048490392850;

#ifdef UTILS_OLD_CAMELCASE
  //! Returns the machine epsilon using camel case style for the specified type.
  //!
  //! \tparam T The type of the machine epsilon. Supported types include float
  //! and double.
  //! \return The machine epsilon constant of type T.
  //! \deprecated
  template <typename T>
  inline T
  machineEps()
  {
    return machine_eps<T>();
  }

  //! Returns the square root of machine epsilon using camel case style for the
  //! specified type.
  //!
  //! \tparam T The type for which to compute the square root of machine
  //! epsilon. Supported types include float and double.
  //! \return The square root of the machine epsilon of type T.
  //! \deprecated
  template <typename T>
  inline T
  sqrtMachineEps()
  {
    return sqrt_machine_eps<T>();
  }

  //! Returns the maximum representable value using camel case style for the
  //! specified type.
  //!
  //! \tparam T The type for which to compute the maximum value. Supported types
  //! include float and double.
  //! \return The maximum representable value of type T.
  //! \deprecated
  template <typename T>
  inline T
  maximumValue()
  {
    return maximum_value<T>();
  }

  //! Returns the minimum representable value using camel case style for the
  //! specified type.
  //!
  //! \tparam T The type for which to compute the minimum value. Supported types
  //! include float and double.
  //! \return The minimum representable value of type T.
  //! \deprecated
  template <typename T>
  inline T
  minimumValue()
  {
    return minimum_value<T>();
  }

  //! Checks if a double value is zero using camel case style.
  //!
  //! \param x The double value to check.
  //! \return True if x is zero, otherwise false.
  //! \deprecated
  static inline bool
  isZero( double x )
  {
    return is_zero( x );
  }

  //! Checks if a float value is zero using camel case style.
  //!
  //! \param x The float value to check.
  //! \return True if x is zero, otherwise false.
  //! \deprecated
  static inline bool
  isZero( float x )
  {
    return is_zero( x );
  }

  //! Checks if a double value is infinite using camel case style.
  //!
  //! \param x The double value to check.
  //! \return True if x is infinite, otherwise false.
  //! \deprecated
  static inline bool
  isInfinite( double x )
  {
    return is_infinite( x );
  }

  //! Checks if a float value is infinite using camel case style.
  //!
  //! \param x The float value to check.
  //! \return True if x is infinite, otherwise false.
  //! \deprecated
  static inline bool
  isInfinite( float x )
  {
    return is_infinite( x );
  }

  //! Checks if a double value is NaN using camel case style.
  //!
  //! \param x The double value to check.
  //! \return True if x is NaN, otherwise false.
  //! \deprecated
  static inline bool
  isNaN( double x )
  {
    return is_NaN( x );
  }

  //! Checks if a float value is NaN using camel case style.
  //!
  //! \param x The float value to check.
  //! \return True if x is NaN, otherwise false.
  //! \deprecated
  static inline bool
  isNaN( float x )
  {
    return is_NaN( x );
  }

  //! Checks if a double value is finite using camel case style.
  //!
  //! \param x The double value to check.
  //! \return True if x is finite, otherwise false.
  //! \deprecated
  static inline bool
  isFinite( double x )
  {
    return is_finite( x );
  }

  //! Checks if a float value is finite using camel case style.
  //!
  //! \param x The float value to check.
  //! \return True if x is finite, otherwise false.
  //! \deprecated
  static inline bool
  isFinite( float x )
  {
    return is_finite( x );
  }

  //! Checks if a double value is a regular (finite) number using camel case
  //! style.
  //!
  //! \param x The double value to check.
  //! \return True if x is a regular number, otherwise false.
  //! \deprecated
  static inline bool
  isRegular( double x )
  {
    return is_finite( x );
  }

  //! Checks if a float value is a regular (finite) number using camel case
  //! style.
  //!
  //! \param x The float value to check.
  //! \return True if x is a regular number, otherwise false.
  //! \deprecated
  static inline bool
  isRegular( float x )
  {
    return is_finite( x );
  }

  //! Checks if a double value is an integer using camel case style.
  //!
  //! \param x The double value to check.
  //! \return True if x is an integer, otherwise false.
  //! \deprecated
  static inline bool
  isInteger( double x )
  {
    return is_integer( x );
  }

  //! Checks if a float value is an integer using camel case style.
  //!
  //! \param x The float value to check.
  //! \return True if x is an integer, otherwise false.
  //! \deprecated
  static inline bool
  isInteger( float x )
  {
    return is_integer( x );
  }

  //! Checks if a double value is unsigned (non-negative integer) using camel
  //! case style.
  //!
  //! \param x The double value to check.
  //! \return True if x is unsigned, otherwise false.
  //! \deprecated
  static inline bool
  isUnsigned( double x )
  {
    return is_unsigned( x );
  }

  //! Checks if a float value is unsigned (non-negative integer) using camel
  //! case style.
  //!
  //! \param x The float value to check.
  //! \return True if x is unsigned, otherwise false.
  //! \deprecated
  static inline bool
  isUnsigned( float x )
  {
    return is_unsigned( x );
  }

  //! Checks if a NaN value is found in an array of doubles using camel case
  //! style.
  //!
  //! \param pv Pointer to the array of double values.
  //! \param DIM The dimension (size) of the array.
  //! \return True if a NaN value is found, otherwise false.
  //! \deprecated
  static inline bool
  foundNaN( double const * pv, int DIM )
  {
    return found_NaN( pv, DIM );
  }

  //! Checks if a NaN value is found in an array of floats using camel case
  //! style.
  //!
  //! \param pv Pointer to the array of float values.
  //! \param DIM The dimension (size) of the array.
  //! \return True if a NaN value is found, otherwise false.
  //! \deprecated
  static inline bool
  foundNaN( float const * pv, int DIM )
  {
    return found_NaN( pv, DIM );
  }

  //! Checks for NaN values in an array of doubles and logs an error if found
  //! using camel case style.
  //!
  //! \param pv Pointer to the array of double values.
  //! \param v_name The name of the variable for logging.
  //! \param DIM The dimension (size) of the array.
  //! \param line The line number where the check is performed.
  //! \param file The name of the file where the check is performed.
  //! \deprecated
  static inline void
  checkNaN( double const * pv, string_view v_name, int DIM, int line, string_view file )
  {
    check_NaN( pv, v_name, DIM, line, file );
  }

  //! Checks for NaN values in an array of floats and logs an error if found
  //! using camel case style.
  //!
  //! \param pv Pointer to the array of float values.
  //! \param v_name The name of the variable for logging.
  //! \param DIM The dimension (size) of the array.
  //! \param line The line number where the check is performed.
  //! \param file The name of the file where the check is performed.
  //! \deprecated
  static inline void
  checkNaN( float const * pv, string_view v_name, int DIM, int line, string_view file )
  {
    check_NaN( pv, v_name, DIM, line, file );
  }

#endif
}  // namespace Utils

//
// eof: Number.hxx
//
