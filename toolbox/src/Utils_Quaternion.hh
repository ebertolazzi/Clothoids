/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2011-2024                                                 |
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
// file: Utils_Quaternion.hh
//

#pragma once
#ifndef UTILS_QUATERNION_HH
#define UTILS_QUATERNION_HH

#include "Utils.hh"

namespace Utils
{

  //!
  //! \brief Quaternion class for 3D rotations and orientations.
  //!
  //! \tparam T The floating-point type used for quaternion components.
  //!
  //! \section quaternion_intro Introduction
  //!
  //! Quaternions are an extension of complex numbers to four dimensions,
  //! represented as \f$ Q = a + b\mathbf{i} + c\mathbf{j} + d\mathbf{k} \f$
  //! where \f$ a, b, c, d \in \mathbb{R} \f$ and \f$ \mathbf{i}, \mathbf{j}, \mathbf{k} \f$
  //! are the fundamental quaternion units satisfying:
  //!
  //! \f{align*}{
  //! \mathbf{i}^2 = \mathbf{j}^2 = \mathbf{k}^2 = \mathbf{i}\mathbf{j}\mathbf{k} = -1
  //! \f}
  //!
  //! This class provides comprehensive quaternion operations with a focus on
  //! 3D rotation representations. Quaternions avoid gimbal lock and provide
  //! efficient interpolation compared to Euler angles.
  //!
  //! \section quaternion_rotation Rotation Representation
  //!
  //! A unit quaternion can represent a rotation of angle \f$ \theta \f$ around
  //! axis \f$ \mathbf{u} = (u_x, u_y, u_z) \f$ as:
  //!
  //! \f[
  //! Q = \cos\left(\frac{\theta}{2}\right) + \sin\left(\frac{\theta}{2}\right)(u_x\mathbf{i} + u_y\mathbf{j} +
  //! u_z\mathbf{k})
  //! \f]
  //!
  //! \section quaternion_operations Operations
  //!
  //! The class supports:
  //! - Basic arithmetic (addition, subtraction, multiplication)
  //! - Conjugation and inversion
  //! - Norm and normalization
  //! - Rotation of vectors
  //! - Conversion to/from axis-angle, rotation matrix, and Euler angles
  //! - Spherical linear interpolation (SLERP)
  //! - Exponential and logarithmic maps
  //!
  //! \section quaternion_references References
  //!
  //! 1. J. B. Kuipers, *Quaternions and Rotation Sequences*, Princeton University Press, 1999.
  //! 2. E. B. Dam, M. Koch, M. Lillholm, *Quaternions, Interpolation and Animation*, Datalogisk Institut, 1998.
  //! 3. Ken Shoemake, *Animating Rotation with Quaternion Curves*, SIGGRAPH 1985.
  //!
  template <typename T> class Quaternion
  {
  public:
    //!
    //! \brief Type definition for the real number type of the quaternion.
    //!
    using real_type = T;

  private:
    real_type m_Q[4];  ///< Quaternion components: [scalar, i, j, k]

  public:
    // ========================================================================
    // CONSTRUCTORS AND BASIC OPERATIONS
    // ========================================================================

    //!
    //! \brief Default constructor (identity quaternion).
    //!
    //! Creates an identity quaternion \f$ Q = 1 + 0\mathbf{i} + 0\mathbf{j} + 0\mathbf{k} \f$
    //! representing no rotation.
    //!
    Quaternion() : m_Q{ 1, 0, 0, 0 } {}

    //!
    //! \brief Constructs a quaternion with specified components.
    //!
    //! \param a The scalar part \f$ a \f$.
    //! \param b The \f$ \mathbf{i} \f$ component \f$ b \f$.
    //! \param c The \f$ \mathbf{j} \f$ component \f$ c \f$.
    //! \param d The \f$ \mathbf{k} \f$ component \f$ d \f$.
    //!
    Quaternion( real_type a, real_type b, real_type c, real_type d ) : m_Q{ a, b, c, d } {}

    //!
    //! \brief Constructs a quaternion from an array.
    //!
    //! \param q Array of 4 elements [scalar, i, j, k].
    //!
    explicit Quaternion( const real_type q[4] ) : m_Q{ q[0], q[1], q[2], q[3] } {}

    //!
    //! \brief Copy constructor.
    //!
    Quaternion( const Quaternion & other ) = default;

    //!
    //! \brief Move constructor.
    //!
    Quaternion( Quaternion && other ) = default;

    //!
    //! \brief Copy assignment operator.
    //!
    Quaternion & operator=( const Quaternion & other ) = default;

    //!
    //! \brief Move assignment operator.
    //!
    Quaternion & operator=( Quaternion && other ) = default;

    // ========================================================================
    // STATIC FACTORY METHODS
    // ========================================================================

    //!
    //! \brief Creates an identity quaternion.
    //!
    //! \return Identity quaternion \f$ Q = 1 + 0\mathbf{i} + 0\mathbf{j} + 0\mathbf{k} \f$.
    //!
    static Quaternion Identity() { return Quaternion( 1, 0, 0, 0 ); }

    //!
    //! \brief Creates a quaternion from axis-angle representation.
    //!
    //! \param axis The rotation axis (must be normalized).
    //! \param angle The rotation angle in radians.
    //! \return The corresponding quaternion.
    //!
    //! \note The axis should be normalized. If not, it will be normalized internally.
    //!
    static Quaternion FromAxisAngle( const real_type axis[3], real_type angle )
    {
      real_type half_angle = angle * real_type( 0.5 );
      real_type sin_half   = std::sin( half_angle );
      real_type cos_half   = std::cos( half_angle );

      // Normalize axis
      real_type norm     = std::sqrt( axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2] );
      real_type inv_norm = norm > real_type( 0 ) ? real_type( 1 ) / norm : real_type( 0 );

      return Quaternion(
        cos_half,
        axis[0] * inv_norm * sin_half,
        axis[1] * inv_norm * sin_half,
        axis[2] * inv_norm * sin_half );
    }

    //!
    //! \brief Creates a quaternion from Euler angles (Z-Y-X convention).
    //!
    //! Uses the aerospace convention: yaw (Z), pitch (Y), roll (X).
    //!
    //! \param yaw Rotation around Z axis (radians).
    //! \param pitch Rotation around Y axis (radians).
    //! \param roll Rotation around X axis (radians).
    //! \return The corresponding quaternion.
    //!
    static Quaternion FromEulerAngles( real_type yaw, real_type pitch, real_type roll )
    {
      real_type cy = std::cos( yaw * real_type( 0.5 ) );
      real_type sy = std::sin( yaw * real_type( 0.5 ) );
      real_type cp = std::cos( pitch * real_type( 0.5 ) );
      real_type sp = std::sin( pitch * real_type( 0.5 ) );
      real_type cr = std::cos( roll * real_type( 0.5 ) );
      real_type sr = std::sin( roll * real_type( 0.5 ) );

      return Quaternion(
        cr * cp * cy + sr * sp * sy,
        sr * cp * cy - cr * sp * sy,
        cr * sp * cy + sr * cp * sy,
        cr * cp * sy - sr * sp * cy );
    }

    //!
    //! \brief Creates a quaternion from a rotation matrix.
    //!
    //! \param R 3x3 rotation matrix (row-major).
    //! \return The corresponding quaternion.
    //!
    //! \note The matrix must be a proper rotation matrix (orthogonal with determinant +1).
    //!
    static Quaternion FromRotationMatrix( const real_type R[3][3] )
    {
      real_type trace = R[0][0] + R[1][1] + R[2][2];
      real_type q[4];

      if ( trace > real_type( 0 ) )
      {
        real_type s = real_type( 0.5 ) / std::sqrt( trace + real_type( 1 ) );
        q[0]        = real_type( 0.25 ) / s;
        q[1]        = ( R[2][1] - R[1][2] ) * s;
        q[2]        = ( R[0][2] - R[2][0] ) * s;
        q[3]        = ( R[1][0] - R[0][1] ) * s;
      }
      else
      {
        if ( R[0][0] > R[1][1] && R[0][0] > R[2][2] )
        {
          real_type s = real_type( 2 ) * std::sqrt( real_type( 1 ) + R[0][0] - R[1][1] - R[2][2] );
          q[0]        = ( R[2][1] - R[1][2] ) / s;
          q[1]        = real_type( 0.25 ) * s;
          q[2]        = ( R[0][1] + R[1][0] ) / s;
          q[3]        = ( R[0][2] + R[2][0] ) / s;
        }
        else if ( R[1][1] > R[2][2] )
        {
          real_type s = real_type( 2 ) * std::sqrt( real_type( 1 ) + R[1][1] - R[0][0] - R[2][2] );
          q[0]        = ( R[0][2] - R[2][0] ) / s;
          q[1]        = ( R[0][1] + R[1][0] ) / s;
          q[2]        = real_type( 0.25 ) * s;
          q[3]        = ( R[1][2] + R[2][1] ) / s;
        }
        else
        {
          real_type s = real_type( 2 ) * std::sqrt( real_type( 1 ) + R[2][2] - R[0][0] - R[1][1] );
          q[0]        = ( R[1][0] - R[0][1] ) / s;
          q[1]        = ( R[0][2] + R[2][0] ) / s;
          q[2]        = ( R[1][2] + R[2][1] ) / s;
          q[3]        = real_type( 0.25 ) * s;
        }
      }

      return Quaternion( q[0], q[1], q[2], q[3] );
    }

    // ========================================================================
    // COMPONENT ACCESS AND MODIFICATION
    // ========================================================================

    //!
    //! \brief Accesses the \f$i\f$-th component of the quaternion.
    //!
    //! \param i Index of the component (0: scalar, 1: i, 2: j, 3: k).
    //! \return The value of the specified component.
    //!
    real_type operator[]( int i ) const { return m_Q[i]; }

    //!
    //! \brief Accesses the \f$i\f$-th component of the quaternion (non-const).
    //!
    //! \param i Index of the component (0: scalar, 1: i, 2: j, 3: k).
    //! \return Reference to the specified component.
    //!
    real_type & operator[]( int i ) { return m_Q[i]; }

    //!
    //! \brief Gets the scalar (real) part of the quaternion.
    //!
    //! \return The scalar part \f$ a \f$.
    //!
    real_type scalar() const { return m_Q[0]; }

    //!
    //! \brief Gets the vector (imaginary) part of the quaternion.
    //!
    //! \param vec Array to store the vector part \f$ [b, c, d] \f$.
    //!
    void vector( real_type vec[3] ) const
    {
      vec[0] = m_Q[1];
      vec[1] = m_Q[2];
      vec[2] = m_Q[3];
    }

    //!
    //! \brief Sets the components of the quaternion.
    //!
    //! \param a The scalar part \f$ a \f$.
    //! \param b The \f$ \mathbf{i} \f$ component \f$ b \f$.
    //! \param c The \f$ \mathbf{j} \f$ component \f$ c \f$.
    //! \param d The \f$ \mathbf{k} \f$ component \f$ d \f$.
    //!
    void set( real_type a, real_type b, real_type c, real_type d )
    {
      m_Q[0] = a;
      m_Q[1] = b;
      m_Q[2] = c;
      m_Q[3] = d;
    }

    //!
    //! \brief Sets the quaternion to identity (no rotation).
    //!
    void set_identity()
    {
      m_Q[0] = 1;
      m_Q[1] = 0;
      m_Q[2] = 0;
      m_Q[3] = 0;
    }

    // ========================================================================
    // QUATERNION OPERATIONS
    // ========================================================================

    //!
    //! \brief Computes the conjugate of the quaternion.
    //!
    //! The conjugate of \f$ Q = a + b\mathbf{i} + c\mathbf{j} + d\mathbf{k} \f$ is:
    //! \f[ \overline{Q} = a - b\mathbf{i} - c\mathbf{j} - d\mathbf{k} \f]
    //!
    //! \return The conjugate quaternion.
    //!
    Quaternion conjugate() const { return Quaternion( m_Q[0], -m_Q[1], -m_Q[2], -m_Q[3] ); }

    //!
    //! \brief Conjugates this quaternion in-place.
    //!
    void conj()
    {
      m_Q[1] = -m_Q[1];
      m_Q[2] = -m_Q[2];
      m_Q[3] = -m_Q[3];
    }

    //!
    //! \brief Computes the inverse of the quaternion.
    //!
    //! The inverse of \f$ Q \f$ is:
    //! \f[ Q^{-1} = \frac{\overline{Q}}{\|Q\|^2} \f]
    //!
    //! \return The inverse quaternion.
    //!
    Quaternion inverse() const
    {
      real_type norm2 = m_Q[0] * m_Q[0] + m_Q[1] * m_Q[1] + m_Q[2] * m_Q[2] + m_Q[3] * m_Q[3];
      if ( norm2 > real_type( 0 ) )
      {
        real_type inv_norm2 = real_type( 1 ) / norm2;
        return Quaternion( m_Q[0] * inv_norm2, -m_Q[1] * inv_norm2, -m_Q[2] * inv_norm2, -m_Q[3] * inv_norm2 );
      }
      return Quaternion();  // Return identity for zero quaternion
    }

    //!
    //! \brief Inverts this quaternion in-place.
    //!
    void invert()
    {
      real_type norm2 = m_Q[0] * m_Q[0] + m_Q[1] * m_Q[1] + m_Q[2] * m_Q[2] + m_Q[3] * m_Q[3];
      if ( norm2 > real_type( 0 ) )
      {
        real_type inv_norm2 = real_type( 1 ) / norm2;
        m_Q[0] *= inv_norm2;
        m_Q[1] *= -inv_norm2;
        m_Q[2] *= -inv_norm2;
        m_Q[3] *= -inv_norm2;
      }
      else
      {
        set_identity();
      }
    }

    //!
    //! \brief Computes the norm (magnitude) of the quaternion.
    //!
    //! \f[ \|Q\| = \sqrt{a^2 + b^2 + c^2 + d^2} \f]
    //!
    //! \return The norm of the quaternion.
    //!
    real_type norm() const
    {
      return std::sqrt( m_Q[0] * m_Q[0] + m_Q[1] * m_Q[1] + m_Q[2] * m_Q[2] + m_Q[3] * m_Q[3] );
    }

    //!
    //! \brief Computes the squared norm of the quaternion.
    //!
    //! \f[ \|Q\|^2 = a^2 + b^2 + c^2 + d^2 \f]
    //!
    //! \return The squared norm of the quaternion.
    //!
    real_type norm_squared() const { return m_Q[0] * m_Q[0] + m_Q[1] * m_Q[1] + m_Q[2] * m_Q[2] + m_Q[3] * m_Q[3]; }

    //!
    //! \brief Normalizes the quaternion to unit length.
    //!
    //! \return Reference to this quaternion after normalization.
    //!
    //! \note If the quaternion has zero norm, it becomes the identity quaternion.
    //!
    Quaternion & normalize()
    {
      real_type n = norm();
      if ( n > real_type( 0 ) )
      {
        real_type inv_n = real_type( 1 ) / n;
        m_Q[0] *= inv_n;
        m_Q[1] *= inv_n;
        m_Q[2] *= inv_n;
        m_Q[3] *= inv_n;
      }
      else
      {
        set_identity();
      }
      return *this;
    }

    //!
    //! \brief Returns a normalized copy of the quaternion.
    //!
    //! \return Normalized quaternion.
    //!
    Quaternion normalized() const
    {
      Quaternion q( *this );
      return q.normalize();
    }

    //!
    //! \brief Checks if the quaternion is normalized (unit length).
    //!
    //! \param tolerance Tolerance for checking unit length.
    //! \return True if the quaternion is normalized within tolerance.
    //!
    bool isNormalized( real_type tolerance = real_type( 1e-6 ) ) const
    {
      real_type diff = std::abs( norm_squared() - real_type( 1 ) );
      return diff <= tolerance;
    }

    //!
    //! \brief Computes the dot product with another quaternion.
    //!
    //! \f[ Q_1 \cdot Q_2 = a_1 a_2 + b_1 b_2 + c_1 c_2 + d_1 d_2 \f]
    //!
    //! \param other The other quaternion.
    //! \return The dot product.
    //!
    real_type dot( const Quaternion & other ) const
    {
      return m_Q[0] * other.m_Q[0] + m_Q[1] * other.m_Q[1] + m_Q[2] * other.m_Q[2] + m_Q[3] * other.m_Q[3];
    }

    // ========================================================================
    // ROTATION OPERATIONS
    // ========================================================================

    //!
    //! \brief Rotates a 3D vector using this quaternion.
    //!
    //! For a unit quaternion \f$ Q \f$, the rotated vector \f$ \mathbf{v}' \f$
    //! is computed as:
    //! \f[ \mathbf{v}' = Q \mathbf{v} Q^{-1} = Q \mathbf{v} \overline{Q} \f]
    //!
    //! \param v The input vector to rotate.
    //! \param w The output rotated vector.
    //!
    void rotate( const real_type v[3], real_type w[3] ) const
    {
      // Extract components for clarity
      real_type a = m_Q[0];
      real_type b = m_Q[1];
      real_type c = m_Q[2];
      real_type d = m_Q[3];

      real_type t2  = a * b;
      real_type t3  = a * c;
      real_type t4  = a * d;
      real_type t5  = -b * b;
      real_type t6  = b * c;
      real_type t7  = b * d;
      real_type t8  = -c * c;
      real_type t9  = c * d;
      real_type t10 = -d * d;

      w[0] = 2 * ( ( t8 + t10 ) * v[0] + ( t6 - t4 ) * v[1] + ( t3 + t7 ) * v[2] ) + v[0];
      w[1] = 2 * ( ( t4 + t6 ) * v[0] + ( t5 + t10 ) * v[1] + ( t9 - t2 ) * v[2] ) + v[1];
      w[2] = 2 * ( ( t7 - t3 ) * v[0] + ( t2 + t9 ) * v[1] + ( t5 + t8 ) * v[2] ) + v[2];
    }

    //!
    //! \brief Rotates a 3D vector and returns the result.
    //!
    //! \param v The input vector to rotate.
    //! \return The rotated vector.
    //!
    std::array<real_type, 3> rotate( const std::array<real_type, 3> & v ) const
    {
      std::array<real_type, 3> w;
      rotate( v.data(), w.data() );
      return w;
    }

    //!
    //! \brief Converts the quaternion to axis-angle representation.
    //!
    //! \param axis Array to store the rotation axis (normalized).
    //! \return The rotation angle in radians.
    //!
    //! \note For identity quaternion, axis is set to [1, 0, 0].
    //!
    real_type to_axis_angle( real_type axis[3] ) const
    {
      Quaternion q     = normalized();
      real_type  angle = real_type( 2 ) * std::acos( q.m_Q[0] );

      real_type sin_half = std::sin( angle * real_type( 0.5 ) );
      if ( sin_half > real_type( 0 ) )
      {
        real_type inv_sin_half = real_type( 1 ) / sin_half;
        axis[0]                = q.m_Q[1] * inv_sin_half;
        axis[1]                = q.m_Q[2] * inv_sin_half;
        axis[2]                = q.m_Q[3] * inv_sin_half;
      }
      else
      {
        axis[0] = 1;
        axis[1] = 0;
        axis[2] = 0;
      }

      return angle;
    }

    //!
    //! \brief Converts the quaternion to a 3x3 rotation matrix.
    //!
    //! \param R 3x3 array to store the rotation matrix (row-major).
    //!
    void to_rotation_matrix( real_type R[3][3] ) const
    {
      Quaternion q = normalized();
      real_type  a = q.m_Q[0];
      real_type  b = q.m_Q[1];
      real_type  c = q.m_Q[2];
      real_type  d = q.m_Q[3];

      real_type aa = a * a, bb = b * b, cc = c * c, dd = d * d;
      real_type ab = a * b, ac = a * c, ad = a * d;
      real_type bc = b * c, bd = b * d, cd = c * d;

      R[0][0] = aa + bb - cc - dd;
      R[0][1] = 2 * ( bc - ad );
      R[0][2] = 2 * ( bd + ac );

      R[1][0] = 2 * ( bc + ad );
      R[1][1] = aa - bb + cc - dd;
      R[1][2] = 2 * ( cd - ab );

      R[2][0] = 2 * ( bd - ac );
      R[2][1] = 2 * ( cd + ab );
      R[2][2] = aa - bb - cc + dd;
    }

    //!
    //! \brief Converts the quaternion to Euler angles (Z-Y-X convention).
    //!
    //! \param yaw Rotation around Z axis (radians).
    //! \param pitch Rotation around Y axis (radians).
    //! \param roll Rotation around X axis (radians).
    //!
    //! \note Uses the aerospace convention: yaw (Z), pitch (Y), roll (X).
    //!
    void to_Euler_angles( real_type & yaw, real_type & pitch, real_type & roll ) const
    {
      Quaternion q = normalized();

      // Roll (x-axis rotation)
      real_type sinr_cosp = real_type( 2 ) * ( q.m_Q[0] * q.m_Q[1] + q.m_Q[2] * q.m_Q[3] );
      real_type cosr_cosp = real_type( 1 ) - real_type( 2 ) * ( q.m_Q[1] * q.m_Q[1] + q.m_Q[2] * q.m_Q[2] );
      roll                = std::atan2( sinr_cosp, cosr_cosp );

      // Pitch (y-axis rotation)
      real_type sinp = real_type( 2 ) * ( q.m_Q[0] * q.m_Q[2] - q.m_Q[3] * q.m_Q[1] );
      if ( std::abs( sinp ) >= real_type( 1 ) )
        pitch = std::copysign( M_PI / real_type( 2 ), sinp );  // Use 90 degrees if out of range
      else
        pitch = std::asin( sinp );

      // Yaw (z-axis rotation)
      real_type siny_cosp = real_type( 2 ) * ( q.m_Q[0] * q.m_Q[3] + q.m_Q[1] * q.m_Q[2] );
      real_type cosy_cosp = real_type( 1 ) - real_type( 2 ) * ( q.m_Q[2] * q.m_Q[2] + q.m_Q[3] * q.m_Q[3] );
      yaw                 = std::atan2( siny_cosp, cosy_cosp );
    }

    // ========================================================================
    // ADVANCED OPERATIONS
    // ========================================================================

    //!
    //! \brief Computes the exponential of the quaternion.
    //!
    //! For a quaternion \f$ Q = a + \mathbf{v} \f$, where \f$ \mathbf{v} \f$ is the vector part:
    //! \f[ \exp(Q) = \exp(a) \left( \cos(\|\mathbf{v}\|) + \frac{\mathbf{v}}{\|\mathbf{v}\|} \sin(\|\mathbf{v}\|)
    //! \right) \f]
    //!
    //! \return The exponential quaternion.
    //!
    Quaternion exp() const
    {
      real_type a      = m_Q[0];
      real_type v_norm = std::sqrt( m_Q[1] * m_Q[1] + m_Q[2] * m_Q[2] + m_Q[3] * m_Q[3] );

      real_type exp_a = std::exp( a );
      real_type cos_v = std::cos( v_norm );

      if ( v_norm > real_type( 0 ) )
      {
        real_type sin_v = std::sin( v_norm );
        real_type scale = exp_a * sin_v / v_norm;
        return Quaternion( exp_a * cos_v, m_Q[1] * scale, m_Q[2] * scale, m_Q[3] * scale );
      }
      else
      {
        return Quaternion( exp_a * cos_v, 0, 0, 0 );
      }
    }

    //!
    //! \brief Computes the natural logarithm of the quaternion.
    //!
    //! For a unit quaternion \f$ Q = \cos(\theta) + \mathbf{u}\sin(\theta) \f$:
    //! \f[ \ln(Q) = \theta \mathbf{u} \f]
    //!
    //! For a general quaternion \f$ Q = r(\cos(\theta) + \mathbf{u}\sin(\theta)) \f$:
    //! \f[ \ln(Q) = \ln(r) + \theta \mathbf{u} \f]
    //!
    //! \return The logarithm quaternion.
    //!
    Quaternion log() const
    {
      real_type a      = m_Q[0];
      real_type v_norm = std::sqrt( m_Q[1] * m_Q[1] + m_Q[2] * m_Q[2] + m_Q[3] * m_Q[3] );
      real_type q_norm = norm();

      if ( q_norm > real_type( 0 ) && v_norm > real_type( 0 ) )
      {
        real_type scale = std::acos( a / q_norm ) / v_norm;
        return Quaternion( std::log( q_norm ), m_Q[1] * scale, m_Q[2] * scale, m_Q[3] * scale );
      }
      else
      {
        // Logarithm of zero or real quaternion
        return Quaternion( std::log( q_norm ), 0, 0, 0 );
      }
    }

    //!
    //! \brief Computes the power of the quaternion.
    //!
    //! \f[ Q^t = \exp(t \ln(Q)) \f]
    //!
    //! \param t The exponent.
    //! \return The quaternion raised to power t.
    //!
    Quaternion pow( real_type t ) const { return ( t * log() ).exp(); }

    //!
    //! \brief Performs spherical linear interpolation (SLERP) between two quaternions.
    //!
    //! \param q1 The starting quaternion.
    //! \param q2 The ending quaternion.
    //! \param t Interpolation parameter in [0, 1].
    //! \param shortestPath If true, takes the shortest path on the hypersphere.
    //! \return The interpolated quaternion.
    //!
    //! \note Both quaternions should be unit quaternions.
    //!
    static Quaternion slerp( const Quaternion & q1, const Quaternion & q2, real_type t, bool shortestPath = true )
    {
      real_type cos_theta = q1.dot( q2 );

      // If shortest path is needed and cos_theta is negative, flip one quaternion
      Quaternion q2_adj = q2;
      if ( shortestPath && cos_theta < real_type( 0 ) )
      {
        cos_theta = -cos_theta;
        q2_adj    = -q2;
      }

      // If the quaternions are very close, use linear interpolation
      if ( cos_theta > real_type( 0.9995 ) )
      {
        Quaternion result = q1 * ( real_type( 1 ) - t ) + q2_adj * t;
        return result.normalized();
      }

      // Spherical interpolation
      real_type theta     = std::acos( cos_theta );
      real_type sin_theta = std::sin( theta );

      real_type w1 = std::sin( ( real_type( 1 ) - t ) * theta ) / sin_theta;
      real_type w2 = std::sin( t * theta ) / sin_theta;

      return q1 * w1 + q2_adj * w2;
    }

    //!
    //! \brief Performs normalized linear interpolation (NLERP) between two quaternions.
    //!
    //! \param q1 The starting quaternion.
    //! \param q2 The ending quaternion.
    //! \param t Interpolation parameter in [0, 1].
    //! \return The interpolated quaternion.
    //!
    //! \note Faster than SLERP but not constant angular velocity.
    //!
    static Quaternion nlerp( const Quaternion & q1, const Quaternion & q2, real_type t )
    {
      Quaternion result = q1 * ( real_type( 1 ) - t ) + q2 * t;
      return result.normalized();
    }

    // ========================================================================
    // ARITHMETIC OPERATORS
    // ========================================================================

    //!
    //! \brief Quaternion addition.
    //!
    Quaternion operator+( const Quaternion & other ) const
    {
      return Quaternion( m_Q[0] + other.m_Q[0], m_Q[1] + other.m_Q[1], m_Q[2] + other.m_Q[2], m_Q[3] + other.m_Q[3] );
    }

    //!
    //! \brief Quaternion subtraction.
    //!
    Quaternion operator-( const Quaternion & other ) const
    {
      return Quaternion( m_Q[0] - other.m_Q[0], m_Q[1] - other.m_Q[1], m_Q[2] - other.m_Q[2], m_Q[3] - other.m_Q[3] );
    }

    //!
    //! \brief Quaternion negation.
    //!
    Quaternion operator-() const { return Quaternion( -m_Q[0], -m_Q[1], -m_Q[2], -m_Q[3] ); }

    //!
    //! \brief Quaternion multiplication.
    //!
    //! Hamilton product:
    //! \f[
    //! \begin{aligned}
    //! Q_1 Q_2 = &(a_1 a_2 - b_1 b_2 - c_1 c_2 - d_1 d_2) \\
    //!          &+ (a_1 b_2 + b_1 a_2 + c_1 d_2 - d_1 c_2)\mathbf{i} \\
    //!          &+ (a_1 c_2 - b_1 d_2 + c_1 a_2 + d_1 b_2)\mathbf{j} \\
    //!          &+ (a_1 d_2 + b_1 c_2 - c_1 b_2 + d_1 a_2)\mathbf{k}
    //! \end{aligned}
    //! \f]
    //!
    Quaternion operator*( const Quaternion & other ) const
    {
      return Quaternion(
        m_Q[0] * other.m_Q[0] - m_Q[1] * other.m_Q[1] - m_Q[2] * other.m_Q[2] - m_Q[3] * other.m_Q[3],
        m_Q[0] * other.m_Q[1] + m_Q[1] * other.m_Q[0] + m_Q[2] * other.m_Q[3] - m_Q[3] * other.m_Q[2],
        m_Q[0] * other.m_Q[2] - m_Q[1] * other.m_Q[3] + m_Q[2] * other.m_Q[0] + m_Q[3] * other.m_Q[1],
        m_Q[0] * other.m_Q[3] + m_Q[1] * other.m_Q[2] - m_Q[2] * other.m_Q[1] + m_Q[3] * other.m_Q[0] );
    }

    //!
    //! \brief Scalar multiplication.
    //!
    Quaternion operator*( real_type scalar ) const
    {
      return Quaternion( m_Q[0] * scalar, m_Q[1] * scalar, m_Q[2] * scalar, m_Q[3] * scalar );
    }

    //!
    //! \brief Scalar division.
    //!
    Quaternion operator/( real_type scalar ) const
    {
      real_type inv_scalar = real_type( 1 ) / scalar;
      return *this * inv_scalar;
    }

    //!
    //! \brief In-place addition.
    //!
    Quaternion & operator+=( const Quaternion & other )
    {
      m_Q[0] += other.m_Q[0];
      m_Q[1] += other.m_Q[1];
      m_Q[2] += other.m_Q[2];
      m_Q[3] += other.m_Q[3];
      return *this;
    }

    //!
    //! \brief In-place subtraction.
    //!
    Quaternion & operator-=( const Quaternion & other )
    {
      m_Q[0] -= other.m_Q[0];
      m_Q[1] -= other.m_Q[1];
      m_Q[2] -= other.m_Q[2];
      m_Q[3] -= other.m_Q[3];
      return *this;
    }

    //!
    //! \brief In-place multiplication.
    //!
    Quaternion & operator*=( const Quaternion & other )
    {
      *this = *this * other;
      return *this;
    }

    //!
    //! \brief In-place scalar multiplication.
    //!
    Quaternion & operator*=( real_type scalar )
    {
      m_Q[0] *= scalar;
      m_Q[1] *= scalar;
      m_Q[2] *= scalar;
      m_Q[3] *= scalar;
      return *this;
    }

    //!
    //! \brief In-place scalar division.
    //!
    Quaternion & operator/=( real_type scalar )
    {
      real_type inv_scalar = real_type( 1 ) / scalar;
      m_Q[0] *= inv_scalar;
      m_Q[1] *= inv_scalar;
      m_Q[2] *= inv_scalar;
      m_Q[3] *= inv_scalar;
      return *this;
    }

    // ========================================================================
    // COMPARISON OPERATORS
    // ========================================================================

    //!
    //! \brief Equality comparison with tolerance.
    //!
    //! \param other The other quaternion.
    //! \param tolerance The tolerance for comparison.
    //! \return True if quaternions are equal within tolerance.
    //!
    bool equals( const Quaternion & other, real_type tolerance = real_type( 1e-6 ) ) const
    {
      return std::abs( m_Q[0] - other.m_Q[0] ) <= tolerance && std::abs( m_Q[1] - other.m_Q[1] ) <= tolerance &&
             std::abs( m_Q[2] - other.m_Q[2] ) <= tolerance && std::abs( m_Q[3] - other.m_Q[3] ) <= tolerance;
    }

    //!
    //! \brief Equality operator.
    //!
    bool operator==( const Quaternion & other ) const
    {
      return m_Q[0] == other.m_Q[0] && m_Q[1] == other.m_Q[1] && m_Q[2] == other.m_Q[2] && m_Q[3] == other.m_Q[3];
    }

    //!
    //! \brief Inequality operator.
    //!
    bool operator!=( const Quaternion & other ) const { return !( *this == other ); }

    // ========================================================================
    // UTILITIES
    // ========================================================================

    //!
    //! \brief Prints the quaternion to an output stream.
    //!
    //! \param os The output stream.
    //!
    void print( ostream_type & os ) const
    {
      os << "[ " << m_Q[0] << ", " << m_Q[1] << "i, " << m_Q[2] << "j, " << m_Q[3] << "k ]";
    }

    //!
    //! \brief Converts the quaternion to a string representation.
    //!
    //! \return String representation of the quaternion.
    //!
    std::string to_string() const
    {
      std::ostringstream oss;
      print( oss );
      return oss.str();
    }
  };

  // ==========================================================================
  // NON-MEMBER FUNCTIONS
  // ==========================================================================

  //!
  //! \brief Scalar multiplication (scalar * quaternion).
  //!
  template <typename T> inline Quaternion<T> operator*( T scalar, const Quaternion<T> & q )
  {
    return q * scalar;
  }

  //!
  //! \brief Output stream operator for quaternions.
  //!
  template <typename T> inline ostream_type & operator<<( ostream_type & os, const Quaternion<T> & q )
  {
    q.print( os );
    return os;
  }

  //!
  //! \brief Computes the angular distance between two unit quaternions.
  //!
  //! \param q1 First unit quaternion.
  //! \param q2 Second unit quaternion.
  //! \return Angular distance in radians.
  //!
  template <typename T> inline T angular_distance( const Quaternion<T> & q1, const Quaternion<T> & q2 )
  {
    T dot_product = q1.dot( q2 );
    dot_product   = std::max( T( -1 ), std::min( T( 1 ), dot_product ) );  // Clamp to [-1, 1]
    return T( 2 ) * std::acos( std::abs( dot_product ) );
  }

  //!
  //! \brief Computes the spherical centroid of a set of quaternions.
  //!
  //! \param quaternions Vector of quaternions.
  //! \param max_iterations Maximum iterations for the averaging algorithm.
  //! \param tolerance Convergence tolerance.
  //! \return The spherical centroid quaternion.
  //!
  template <typename T> inline Quaternion<T> spherical_centroid(
    const std::vector<Quaternion<T>> & quaternions,
    int                                max_iterations = 100,
    T                                  tolerance      = T( 1e-6 ) )
  {
    if ( quaternions.empty() ) return Quaternion<T>::Identity();

    // Start with the first quaternion
    Quaternion<T> centroid = quaternions[0].normalized();

    for ( int iter = 0; iter < max_iterations; ++iter )
    {
      // Compute tangent vectors from centroid to each quaternion
      Quaternion<T> tangent_sum( 0, 0, 0, 0 );

      for ( const auto & q : quaternions )
      {
        // Compute tangent vector = log(centroid⁻¹ * q)
        Quaternion<T> delta = centroid.conjugate() * q;
        tangent_sum += delta.log();
      }

      // Average tangent vector
      tangent_sum /= T( quaternions.size() );

      // Update centroid: centroid = centroid * exp(tangent_avg)
      centroid = centroid * tangent_sum.exp();
      centroid.normalize();

      // Check convergence
      if ( tangent_sum.norm() < tolerance ) break;
    }

    return centroid;
  }

}  // namespace Utils

#endif  // UTILS_QUATERNION_HH

//
// eof: Utils_Quaternion.hh
//
