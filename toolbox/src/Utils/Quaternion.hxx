/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2011                                                      |
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
/// file: Quaternion.hxx
///

namespace Utils {

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  using std::ostream;
  #endif

  //!
  //! **Quaternion class**
  //!
  //! A quaternion is a quadruplet (A,B,C,D) of real numbers,
  //! which may be written as
  //!
  //! \f[ Q = A + B \mathbf{i} + C\mathbf{j} + D\mathbf{k} \f]
  //!
  //!
  //! **Reference:**
  //!
  //! - | James Foley, Andries van Dam, Steven Feiner, John Hughes,
  //!   | Computer Graphics, Principles and Practice, Second Edition,
  //!   | Addison Wesley, 1990.
  //!
  template <typename T>
  class Quaternion {
  public:
    //!
    //! Type quaternion of coordinate.
    //!
    using real_type = T;

  private:
    real_type m_Q[4];

  public:

    //!
    //! Build null quaternion.
    //!
    Quaternion() {
      m_Q[0] = m_Q[1] = m_Q[2] = m_Q[3] = 0;
    }

    //!
    //! Build quaternion \f$ A + B \mathbf{i} + C \mathbf{j} + D \mathbf{k} \f$.
    //!
    Quaternion(
      real_type A,
      real_type B,
      real_type C,
      real_type D
    ) {
      m_Q[0] = A;
      m_Q[1] = B;
      m_Q[2] = C;
      m_Q[3] = D;
    }

    //!
    //! Build quaternion \f$ A + B \mathbf{i} + C \mathbf{j} + D \mathbf{k} \f$.
    //!
    void
    setup(
      real_type A,
      real_type B,
      real_type C,
      real_type D
    ) {
      m_Q[0] = A;
      m_Q[1] = B;
      m_Q[2] = C;
      m_Q[3] = D;
    }

    //!
    //! Print a quaternion to stream.
    //!
    void
    print( ostream_type & os) const {
      os << "[ "
         << m_Q[0] << ", "
         << m_Q[1] << "i, "
         << m_Q[2] << "j, "
         << m_Q[3] << "k ]";
    }

    //!
    //! Return `i`-th coordinate of the quaternion
    //!
    real_type
    operator [] (int i) const { return m_Q[i]; }

    //!
    //! Conjugates a quaternion
    //!
    //! \f[ Q = A + B \mathbf{i} + C \mathbf{j} + D \mathbf{k} \f]
    //!
    //! The conjugate of \f$ Q \f$ is
    //!
    //! \f[ \overline{Q} = A - B \mathbf{i} - C \mathbf{j} - D \mathbf{k} \f]
    //!
    void
    conj() { m_Q[1] = -m_Q[1]; m_Q[2] = -m_Q[2]; m_Q[3] = -m_Q[3]; }

    //!
    //! Invert a quaternion \f$ Q = A + B \mathbf{i} + C \mathbf{j} + D \mathbf{k} \f$.
    //! The inverse of of \f$ Q \f$ is
    //! \f[
    //!   Q^{-1} = \dfrac{ A - B \mathbf{i} - C \mathbf{j} - D \mathbf{k}}
    //!                  { A^2 + B^2 + C^2 + D^2 }
    //! \f]
    //!
    void
    invert() {
      real_type bf = 1/(m_Q[0]*m_Q[0] + m_Q[1]*m_Q[1] + m_Q[2]*m_Q[2] + m_Q[3]*m_Q[3]);
      m_Q[0] *= bf; bf = -bf; m_Q[1] *= bf; m_Q[2] *= bf; m_Q[3] *= bf;
    }

    //!
    //! Computes the norm of a quaternion.
    //! The norm of \f$ Q \f$ is \f$ \sqrt{A^2+B^2+C^2+D^2} \f$.
    //!
    real_type
    norm() const {
      return sqrt(m_Q[0]*m_Q[0] + m_Q[1]*m_Q[1] + m_Q[2]*m_Q[2] + m_Q[3]*m_Q[3]);
    }

    //!
    //! Applies a quaternion rotation to a vector in 3D.
    //!
    //! If \f$ Q \f$ is a unit quaternion that encodes a rotation of ANGLE
    //! radians about the vector AXIS, then for an arbitrary real
    //! vector \f$ V \f$, the result \f$ W \f$ of the rotation on \f$ V \f$
    //! can be written as:
    //!
    //! \f$ W = Q * V * Conj(Q) \f$
    //!
    void
    rotate( real_type const v[3], real_type w[3] ) const {
      w[0] = ( m_Q[0] * m_Q[0] + m_Q[1] * m_Q[1] ) * v[0]
           + ( m_Q[1] * m_Q[2] - m_Q[0] * m_Q[3] ) * v[1]
           + ( m_Q[1] * m_Q[3] + m_Q[0] * m_Q[2] ) * v[2];

      w[1] = ( m_Q[1] * m_Q[2] + m_Q[0] * m_Q[3] ) * v[0]
           + ( m_Q[0] * m_Q[0] + m_Q[2] * m_Q[2] ) * v[1]
           + ( m_Q[2] * m_Q[3] - m_Q[0] * m_Q[1] ) * v[2];

      w[2] = ( m_Q[1] * m_Q[3] - m_Q[0] * m_Q[2] ) * v[0]
           + ( m_Q[2] * m_Q[3] + m_Q[0] * m_Q[1] ) * v[1]
           + ( m_Q[0] * m_Q[0] + m_Q[3] * m_Q[3] ) * v[2];

      w[0] = 2*w[0] - v[0];
      w[1] = 2*w[1] - v[1];
      w[2] = 2*w[2] - v[2];
    }

    //!
    //! Converts a rotation from quaternion to axis format in 3D.
    //!
    //!  A rotation quaternion Q has the form:
    //!
    //! \f[ Q = A + B \mathbf{i} + C \mathbf{j} + D \mathbf{k} \f]
    //!
    //! where A, B, C and D are real numbers, and \f$ \mathbf{i} \f$,
    //! \f$ \mathbf{j} \f$, and \f$ \mathbf{k} \f$ are to be regarded
    //! as symbolic constant basis vectors, similar to the role of the
    //! "\f$ \mathbf{i} \f$" in the representation of imaginary numbers.
    //!
    //! A is the cosine of half of the angle of rotation. (B,C,D) is a
    //! vector pointing in the direction of the axis of rotation.
    //! Rotation multiplication and inversion can be carried out using
    //! this format and the usual rules for quaternion multiplication
    //! and inversion.
    //!
    real_type
    to_axis( real_type axis[3] ) const {
      real_type sin_phi = sqrt( m_Q[1]*m_Q[1] + m_Q[2]*m_Q[2] + m_Q[3]*m_Q[3] );
      real_type cos_phi = m_Q[0];
      real_type angle   = 2 * atan2( sin_phi, cos_phi );
      if ( sin_phi == 0 ) {
        axis[0] = 1; axis[1] = axis[2] = 0;
      } else {
        axis[0] = m_Q[1] / sin_phi;
        axis[1] = m_Q[2] / sin_phi;
        axis[2] = m_Q[3] / sin_phi;
      }
      return angle;
    }

    //!
    //! Converts a rotation from quaternion to rotation matrix.
    //!
    void
    to_matrix( real_type mat[3][3] ) const {
      real_type axis[3];
      real_type angle = to_axis( axis );
      real_type ca    = cos( angle );
      real_type sa    = sin( angle );

      mat[0][0] =              axis[0] * axis[0] + ca * ( 1 - axis[0] * axis[0] );
      mat[1][0] = ( 1 - ca ) * axis[0] * axis[1] - sa * axis[2];
      mat[2][0] = ( 1 - ca ) * axis[0] * axis[2] + sa * axis[1];

      mat[0][1] = ( 1 - ca ) * axis[1] * axis[0] + sa * axis[2];
      mat[1][1] =              axis[1] * axis[1] + ca * ( 1 - axis[1] * axis[1] );
      mat[2][1] = ( 1 - ca ) * axis[1] * axis[2] - sa * axis[0];

      mat[0][2] = ( 1 - ca ) * axis[2] * axis[0] - sa * axis[1];
      mat[1][2] = ( 1 - ca ) * axis[2] * axis[1] + sa * axis[0];
      mat[2][2] =              axis[2] * axis[2] + ca * ( 1 - axis[2] * axis[2] );
    }

  };

  //!
  //! Multiplies two quaternions.
  //!
  //! To multiply two quaternions, use the relationships:
  //!
  //! \f[ \mathbf{i} * \mathbf{j} = -\mathbf{j} * \mathbf{i} = \mathbf{k} \f]
  //! \f[ \mathbf{j} * \mathbf{k} = -\mathbf{k} * \mathbf{j} = \mathbf{i} \f]
  //! \f[ \mathbf{k} * \mathbf{i} = -\mathbf{i} * \mathbf{k} = \mathbf{j} \f]
  //! \f[ \mathbf{i} * \mathbf{i} =  \mathbf{j} * \mathbf{j} = \mathbf{k} * \mathbf{k} = -1 \f]
  //!
  template <typename T>
  inline
  Quaternion<T>
  operator * ( Quaternion<T> const & a, Quaternion<T> const & b ) {
    return Quaternion<T>( a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
                          a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
                          a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
                          a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0] );
  }

  template <typename T>
  inline
  ostream_type& operator << ( ostream_type & os,  Quaternion<T> const & Q ) {
    Q.print(os);
    return os;
  }
}

///
/// eof: Quaternion.hxx
///
