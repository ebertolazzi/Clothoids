/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2022                                                      |
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
// file: Utils_Poly.hh
// Header-only polynomial and Sturm sequence library
//
// References:
// [1] Horner, W. G. (1819). "A new method of solving numerical equations
//     of all orders, by continuous approximation".
// [2] Sturm, C. (1829). "Mémoire sur la résolution des équations numériques".
//     Bulletin des Sciences de Férussac.
// [3] Uspensky, J. V. (1948). "Theory of Equations".
// [4] Press, W. H., et al. (2007). "Numerical Recipes: The Art of Scientific
//     Computing", 3rd ed., Cambridge University Press.
//

#pragma once

#ifndef UTILS_POLY_dot_HH
#define UTILS_POLY_dot_HH

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cassert>
#include <limits>

#include "Utils_eigen.hh"
#include "Utils_AlgoBracket.hh"
#include "Utils_fmt.hh"

#if defined( __llvm__ ) || defined( __clang__ )
#pragma clang diagnostic ignored "-Wdeprecated-copy-with-dtor"
#endif

namespace Utils
{

  using std::string;
  using std::vector;

  /**
   * @brief Computes machine epsilon for the given floating-point type.
   *
   * @tparam T Floating-point type
   * @return T Machine epsilon
   */
  template <typename T> inline T machine_eps()
  {
    return std::numeric_limits<T>::epsilon();
  }

  //!
  //! \class Poly
  //! \brief Specializes the `Eigen::Matrix` class to represent and manipulate
  //! polynomials.
  //!
  //! This class provides a convenient way to handle polynomials using the
  //! `Eigen::Matrix` class, where a polynomial is represented as a vector of
  //! coefficients in increasing order (constant term first).
  //!
  //! Given a coefficient vector:
  //!
  //! \f[
  //!   \boldsymbol{v} = ( v_0, v_1, v_2, \ldots, v_m )
  //! \f]
  //!
  //! It corresponds to the polynomial:
  //!
  //! \f[
  //!   p(x) = v_0 + v_1\,x + v_2\,x^2 + \cdots + v_m\,x^m
  //! \f]
  //!
  //! The class overloads common operators (`+`, `-`, `*`) to support polynomial
  //! arithmetic, making it simple to add, subtract, and multiply polynomials.
  //!
  //! **Features**
  //!
  //! - **Operator Overloading**: Supports the `+`, `-`, and `*` operators for
  //! polynomial addition, subtraction, and multiplication, respectively.
  //! - **Flexible Degree Management**: Allows setting and adjusting the polynomial degree.
  //! - **Integration with Eigen**: Leverages the `Eigen::Matrix` class for efficient vector
  //!   and matrix operations.
  //! - **Numerical Methods**: Includes Horner's method for polynomial evaluation
  //!   (see [1]) and Sturm sequence for root isolation (see [2]).
  //!
  //! **Usage**
  //!
  //! \code{cpp}
  //! #include "Utils_Poly.hh"
  //!
  //! Poly<double> P, Q, R;
  //!
  //! // Initialize P(x) = 1 + x + 2*x^2
  //! P.set_degree(2);
  //! P << 1, 1, 2;
  //!
  //! // Initialize Q(x) = 1 - 3*x^4
  //! Q.set_degree(4);
  //! Q << 1, 0, 0, 0, -3;
  //!
  //! // Compute R(x) = P(x) + 3 * Q(x)
  //! R = P + 3 * Q;
  //!
  //! // Compute R(x) = P(x) * Q(x)
  //! R = P * Q;
  //! \endcode
  //!
  //! @tparam Real Floating-point type for coefficients (float, double, etc.)
  //!
  template <typename Real> class Poly : public Eigen::Matrix<Real, Eigen::Dynamic, 1>
  {
  public:
    using Integer = int;
    using dvec_t  = Eigen::Matrix<Real, Eigen::Dynamic, 1>;
    using Poly_t  = Poly<Real>;

  private:
    Integer  m_order;
    dvec_t & to_eigen() { return *static_cast<dvec_t *>( this ); }

  public:
    Poly() : m_order( 0 ) {}

    //!
    //! Initializes the polynomial with a specified maximum degree.
    //!
    //! \param order The highest degree of the polynomial that can be stored.
    //!
    //! This defines the size of the coefficient vector.
    //!
    explicit Poly( int order ) : m_order( order )
    {
      this->resize( order );
      this->setZero();
    }

    //!
    //! Copy constructor for the polynomial.
    //! Initializes the polynomial by creating a copy of the given polynomial
    //! \f$c(x)\f$.
    //!
    //! \param c The polynomial to be copied.
    //!
    Poly( Poly_t const & c ) : m_order( c.m_order ) { *this = c; }

    //!
    //! Constructs a polynomial by copying the given vector of
    //! coefficients \f$c(x)\f$.
    //!
    //! \param c A vector of coefficients representing the polynomial,
    //!        where each element corresponds to the coefficient of
    //!        increasing powers of \f$x\f$ (i.e., `c[i]` is the
    //!        coefficient of \f$x^i\f$).
    //!
    //! The size of the vector determines the polynomial's degree.
    //! The copy operation explicitly invokes the base class
    //! assignment operator to avoid recursion.
    //!
    explicit Poly( dvec_t const & c )
    {
      this->resize( c.size() );
      // per evitare la ricorsione devo chiamare esplicitamente
      // l'operatore = della classe di base.
      this->dvec_t::operator=( c );
      m_order = Integer( c.size() );
    }

    //!
    //! Returns the polynomial as an `Eigen` vector.
    //! This function provides access to the underlying
    //! `Eigen` vector representation of the polynomial's coefficients.
    //!
    //! \return A constant reference to the `Eigen` vector that stores
    //!         the polynomial's coefficients.
    //!          The vector elements correspond to the coefficients of
    //!          increasing powers of \f$x\f$ (i.e., the element at index
    //!          \f$i\f$ is the coefficient of \f$x^i\f$).
    //!
    dvec_t const & to_eigen() const { return *static_cast<dvec_t const *>( this ); }

    //!
    //! Sets the order (degree+1) of the polynomial and resizes the coefficient
    //! vector accordingly.
    //!
    //! This function adjusts the polynomial to the specified degree by resizing
    //! the underlying coefficient vector to match the new order. All
    //! coefficients are initialized to zero.
    //!
    //! \param order The new order of the polynomial.
    //!        The polynomial will have \f$order\f$ terms
    //!        (from the constant term up to the \f$x^{\text{order-1}}\f$ term).
    //!
    //! **Example**
    //!
    //! \code{cpp}
    //! Poly<double> P;
    //!
    //! // Set the polynomial degree to 3 (P(x) = 0 + 0*x + 0*x^2 + 0*x^3)
    //! P.set_order(4);
    //!
    //! // After setting the degree, we can manually set the coefficients
    //! P << 1, -2, 3, 4;  // P(x) = 1 - 2*x + 3*x^2 + 4*x^3
    //! \endcode
    //!
    void set_order( Integer order )
    {
      m_order = order;
      this->resize( order );
      this->setZero();
    }

    //!
    //! Sets the degree of the polynomial and resizes the coefficient vector
    //! accordingly.
    //!
    //! This function adjusts the polynomial to the specified degree by resizing
    //! the underlying coefficient vector to have \f$degree + 1\f$ terms. All
    //! coefficients are initialized to zero.
    //!
    //! \param degree The highest exponent of \f$x\f$ in the polynomial.
    //!               The polynomial will have \f$degree + 1\f$ terms
    //!               (from the constant term up to the \f$x^{\text{degree}}\f$
    //!               term).
    //!
    //! **Example**
    //!
    //! \code{cpp}
    //! Poly<double> P;
    //!
    //! // Set the polynomial degree to 3 (P(x) = 0 + 0*x + 0*x^2 + 0*x^3)
    //! P.set_degree(3);
    //!
    //! // After setting the degree, we can manually assign the coefficients
    //! P << 1, -2, 3, 4;  // P(x) = 1 - 2*x + 3*x^2 + 4*x^3
    //! \endcode
    //!
    void set_degree( Integer degree )
    {
      m_order = degree + 1;
      this->resize( m_order );
      this->setZero();
    }

    //!
    //! Initializes the polynomial as a scalar and returns a
    //! reference to the internal polynomial.
    //!
    //! This function sets the polynomial to a constant value,
    //! effectively representing the polynomial \f$p(x) = a\f$,
    //! where \f$a\f$ is the provided scalar value.
    //!
    //! \param a The scalar value to initialize the polynomial.
    //!
    //! \return A reference to the internal polynomial,
    //!         allowing for method chaining.
    //!
    //! **Example**
    //!
    //! \code{cpp}
    //! Poly<double> P;
    //!
    //! // Set the polynomial to a scalar value of 5 (P(x) = 5)
    //! P.set_scalar(5);
    //! \endcode
    //!
    Poly_t & set_scalar( Real a )
    {
      this->resize( 1 );
      this->coeffRef( 0 ) = a;
      m_order             = 1;
      return *this;
    }

    //!
    //! Initializes the polynomial as \f$x + a\f$ and returns a
    //! reference to the internal polynomial.
    //!
    //! This function sets the polynomial to a linear form,
    //! representing the polynomial
    //!
    //! \f[ p(x) = x + a \f]
    //!
    //! where \f$a\f$ is the provided constant.
    //!
    //! \param a The constant value to be added to the polynomial.
    //!
    //! \return A reference to the internal polynomial, allowing for method
    //! chaining.
    //!
    //! **Example**
    //!
    //! \code{cpp}
    //! Poly<double> P;
    //!
    //! // Set the polynomial to a monomial of the form x + 3 (P(x) = x + 3)
    //! P.set_monomial(3);
    //! \endcode
    //!
    Poly_t & set_monomial( Real a )
    {
      this->resize( 2 );
      this->coeffRef( 0 ) = a;
      this->coeffRef( 1 ) = 1;
      m_order             = 2;
      return *this;
    }

    //!
    //! Returns a copy of the coefficients of the polynomial as an `Eigen`
    //! vector.
    //!
    //! This function provides a copy of the underlying `Eigen` vector that
    //! contains the polynomial's coefficients. The coefficients correspond to
    //! increasing powers of \f$x\f$, with the element at index \f$i\f$ being
    //! the coefficient of \f$x^i\f$.
    //!
    //! \return A copy of the `Eigen` vector containing the polynomial
    //! coefficients.
    //!
    //! **Example**
    //!
    //! \code{cpp}
    //! Poly<double> P;
    //! P.set_monomial(3); // P(x) = x + 3
    //! dvec_t coefficients = P.coeffs(); // coefficients will contain [3, 1]
    //! \endcode
    //!
    dvec_t coeffs() const { return to_eigen(); }

    //!
    //! Returns the degree of the polynomial.
    //!
    //! This function calculates and returns the degree of the polynomial, which
    //! is defined as the highest exponent of \f$x\f$ in the polynomial. The
    //! degree is equal to
    //! \f$ m_order - 1 \f$, reflecting the number of terms in the polynomial
    //! minus one.
    //!
    //! \return The degree of the polynomial as an integer value.
    //!
    //! **Example**
    //!
    //! \code{cpp}
    //! Poly<double> P;
    //! P.set_monomial(3); // P(x) = x + 3
    //! Integer deg = P.degree(); // deg will be 1
    //! \endcode
    //!
    Integer degree() const { return m_order - 1; }

    //!
    //! Returns the order of the polynomial.
    //!
    //! This function returns the order (size) of the polynomial, which is
    //! defined as the total number of coefficients stored in the polynomial,
    //! including the constant term. The order is equal to \f$ m_order \f$.
    //!
    //! \return The order of the polynomial as an integer value, representing
    //! the total
    //!         number of coefficients.
    //!
    //! **Example**
    //!
    //! \code{cpp}
    //! Poly<double> P;
    //! P.set_monomial(3); // P(x) = x + 3
    //! Integer ord = P.order(); // ord will be 2
    //! \endcode
    //!
    Integer order() const { return m_order; }

    //!
    //! Converts the polynomial to a human-readable string representation.
    //!
    //! This method generates a string that represents the polynomial
    //! in the standard mathematical notation.
    //! The output format will display the polynomial's terms with
    //! appropriate signs and powers of \f$x\f$.
    //!
    //! - If the polynomial is empty (order <= 0), it returns "EMPTY!".
    //!
    //! - If the polynomial consists of a single term (order == 1),
    //!   it returns the coefficient of that term.
    //!
    //! - If all coefficients are zero, it returns "0".
    //!
    //! The polynomial is represented in the form:
    //!
    //! \f[ p(x) = c_0 + c_1 x + c_2 x^2 + \ldots + c_n x^n \f]
    //!
    //! where \f$ c_i \f$ are the coefficients of the polynomial.
    //! The function handles formatting for positive and negative
    //! coefficients, omitting terms with a coefficient of zero,
    //! and does not display a coefficient of one when writing
    //! the variable term.
    //!
    //! \return A string representation of the polynomial.
    //!
    //! **Example**
    //!
    //! \code{cpp}
    //! Poly<double> P;
    //! P.set_monomial(3); // P(x) = x + 3
    //! string str = P.to_string(); // str will be "x + 3"
    //! \endcode
    //!
    string to_string() const;

    //!
    //! Evaluates the polynomial at a given point using Horner's method [1].
    //!
    //! This method computes the value of the polynomial \f$p(x)\f$ at
    //! a specified point \f$ x \f$ using Horner's method, which is an
    //! efficient algorithm for polynomial evaluation.
    //! The polynomial is represented in the form:
    //!
    //! \f[ p(x) = c_0 + c_1 x + c_2 x^2 + \ldots + c_n x^n \f]
    //!
    //! where \f$ c_i \f$ are the coefficients of the polynomial.
    //! Horner's method reduces the number of multiplications required
    //! to evaluate the polynomial, providing a more efficient
    //! computation, especially for polynomials of high degree.
    //!
    //! \param x The point at which the polynomial is to be evaluated.
    //!
    //! \return The computed value of the polynomial at \f$ x \f$.
    //!
    //! **Example**
    //!
    //! \code{cpp}
    //! Poly<double> P;
    //! P.set_monomial(2); // P(x) = x + 2
    //! Real value = P.eval(3); // value will be 5 = 3 + 2
    //! \endcode
    //!
    Real eval( Real x ) const
    {
      // Calcolo il polinomio usando il metodo di Horner
      Integer n{ m_order - 1 };
      Real    res{ this->coeff( n ) };
      while ( n-- > 0 ) res = res * x + this->coeff( n );
      return res;
    }

    //!
    //! Evaluates the derivative of the polynomial at a
    //! given point using Horner's method.
    //!
    //! This method computes the value of the derivative \f$ p'(x) \f$
    //! of the polynomial \f$ p(x) \f$ at a specified point \f$ x \f$.
    //! The derivative is calculated using Horner's method, which
    //! is an efficient technique for evaluating polynomials and
    //! their derivatives.
    //!
    //! The polynomial is represented in the form:
    //!
    //! \f[ p(x) = c_0 + c_1 x + c_2 x^2 + \ldots + c_n x^n \f]
    //!
    //! The derivative is given by:
    //!
    //! \f[ p'(x) = c_1 + 2c_2 x + 3c_3 x^2 + \ldots + n c_n x^{n-1} \f]
    //!
    //! where \f$ c_i \f$ are the coefficients of the polynomial.
    //! This implementation efficiently computes the derivative
    //! without explicitly forming it, thus optimizing performance.
    //!
    //! \param x The point at which the derivative of the polynomial is to be
    //! evaluated.
    //!
    //! \return The computed value of the derivative at \f$ x \f$.
    //!
    //! **Example**
    //!
    //! \code{cpp}
    //! Poly<double> P;
    //! P.set_monomial(2); // P(x) = x + 2
    //! Real derivativeValue = P.eval_D(3); // derivativeValue will be 1 (the
    //! derivative is constant)
    //! \endcode
    //!
    Real eval_D( Real x ) const
    {
      // Calcolo il polinomio usando il metodo di Horner
      Integer n{ m_order - 1 };
      Real    Dp{ this->coeff( n ) * n };
      while ( --n > 0 ) Dp = Dp * x + this->coeff( n ) * n;
      return Dp;
    }

    //!
    //! Evaluates the polynomial and its derivative at a given
    //! point using Horner's method.
    //!
    //! This method computes both the value of the polynomial
    //! \f$ p(x) \f$ and its derivative \f$ p'(x) \f$ at a
    //! specified point \f$ x \f$.
    //! The calculations are performed using Horner's
    //! method, which is an efficient algorithm for polynomial
    //! evaluation and its derivatives.
    //!
    //! The polynomial is represented in the form:
    //!
    //! \f[ p(x) = c_0 + c_1 x + c_2 x^2 + \ldots + c_n x^n \f]
    //!
    //! The derivative is represented as:
    //!
    //! \f[ p'(x) = c_1 + 2c_2 x + 3c_3 x^2 + \ldots + n c_n x^{n-1} \f]
    //!
    //! The results are stored in the provided reference parameters, where:
    //! - \f$ p \f$ will contain the value of the polynomial at \f$ x \f$.
    //! - \f$ Dp \f$ will contain the value of the derivative at \f$ x \f$.
    //!
    //! \param x The point at which the polynomial and its derivative are to be
    //! evaluated.
    //!
    //! \param p Reference to a variable where the computed value of the
    //! polynomial
    //!          at \f$ x \f$ will be stored. This variable is modified to hold
    //!          the result of the evaluation.
    //!
    //! \param Dp Reference to a variable where the computed value of the
    //! derivative
    //!           at \f$ x \f$ will be stored. This variable is modified to hold
    //!           the result of the derivative evaluation.
    //!
    //! **Example**
    //!
    //! \code{cpp}
    //! Poly<double> P;
    //! P.set_monomial(2); // P(x) = x + 2
    //! Real value, derivative;
    //! P.eval(3, value, derivative); // value will be 5 and derivative will be
    //! 1
    //! \endcode
    //!
    void eval( Real x, Real & p, Real & Dp ) const
    {
      // Calcolo il polinomio usando il metodo di Horner
      Integer n{ m_order - 1 };
      p  = this->coeff( n );
      Dp = this->coeff( n ) * n;
      while ( --n > 0 )
      {
        p  = p * x + this->coeff( n );
        Dp = Dp * x + this->coeff( n ) * n;
      }
      p = p * x + this->coeff( 0 );
    }

    //!
    //! Returns the leading coefficient of the polynomial.
    //!
    //! The leading coefficient is the coefficient of the term
    //! with the highest degree in the polynomial.
    //! This method retrieves the coefficient corresponding
    //! to the term  \f$ c_n x^n \f$, where  \f$ c_n \f$
    //! is the coefficient stored at the last index of the coefficient vector.
    //!
    //! \return The leading coefficient of the polynomial.
    //!
    //! **Example**
    //!
    //! \code{cpp}
    //! Poly<double> P;
    //! P.set_degree(3); // P(x) = c_0 + c_1 x + c_2 x^2 + c_3 x^3
    //! P << 1, 2, 3, 4; // Coefficients: c_0 = 1, c_1 = 2, c_2 = 3, c_3 = 4
    //! Real leadCoeff = P.leading_coeff(); // leadCoeff will be 4
    //! \endcode
    //!
    Real leading_coeff() const { return this->coeff( m_order - 1 ); }

    //!
    //! Computes the derivative of the polynomial and stores
    //! the result in the provided polynomial object.
    //!
    //! This method calculates the derivative of the polynomial
    //! using the power rule. The derivative of a polynomial
    //! \f$ p(x) = c_0 + c_1 x + c_2 x^2 + \ldots + c_n x^n \f$ is given by:
    //!
    //! \f[ p'(x) = c_1 + 2c_2 x + 3c_3 x^2 + \ldots + n c_n x^{n-1} \f]
    //!
    //! The result is stored in the polynomial object passed as an argument.
    //!
    //! \param result The polynomial object where the derivative will be stored.
    //!
    //! **Example**
    //!
    //! \code{cpp}
    //! Poly<double> P;
    //! P.set_degree(3); // P(x) = c_0 + c_1 x + c_2 x^2 + c_3 x^3
    //! P << 1, 2, 3, 4; // Coefficients: c_0 = 1, c_1 = 2, c_2 = 3, c_3 = 4
    //! Poly<double> derivativePoly;
    //! P.derivative(derivativePoly); // derivativePoly now contains the
    //! derivative
    //! \endcode
    //!
    void derivative( Poly & result ) const
    {
      result.resize( m_order - 1 );  // nuovo polinomio contenente il risultato
      for ( Integer i = 1; i < m_order; ++i ) result.coeffRef( i - 1 ) = i * this->coeff( i );
      result.m_order = m_order - 1;
    }

    //!
    //! Computes the integral of the polynomial and stores the
    //! result in the provided polynomial object.
    //!
    //! This method calculates the indefinite integral of the
    //! polynomial using the power rule.
    //! The integral of a polynomial \f$p(x) = c_0 + c_1 x + c_2 x^2 + \ldots +
    //! c_n x^n \f$ is given by:
    //!
    //! \f[ \int p(x) \,dx = c_0 x + \frac{c_1}{2} x^2 + \frac{c_2}{3} x^3 +
    //! \ldots + \frac{c_n}{n+1} x^{n+1} + C \f]
    //!
    //! The constant of integration  \f$ C \f$ is not included in the result.
    //! The result is stored in the polynomial object passed as an argument.
    //!
    //! \param result The polynomial object where the integral will be stored.
    //!
    //! **Example**
    //!
    //! \code{cpp}
    //! Poly<double> P;
    //! P.set_degree(3); // P(x) = c_0 + c_1 x + c_2 x^2 + c_3 x^3
    //! P << 1, 2, 3, 4; // Coefficients: c_0 = 1, c_1 = 2, c_2 = 3, c_3 = 4
    //! Poly<double> integralPoly;
    //! P.integral(integralPoly); // integralPoly now contains the integral
    //! \endcode
    //!
    void integral( Poly & result ) const
    {
      result.resize( m_order + 1 );  // nuovo polinomio contenente il risultato
      result.coeffRef( 0 ) = 0;
      for ( Integer i = 1; i <= m_order; ++i ) result.coeffRef( i ) = this->coeff( i - 1 ) / i;
      result.m_order = m_order + 1;
    }

    //!
    //! Normalizes the polynomial \f$ p(x) = \sum_{i=0}^n a_i x^i \f$
    //! by scaling its coefficients such that the maximum absolute value
    //! of the coefficients is equal to 1.
    //! The scaling factor \f$ S \f$ is calculated to achieve this
    //! normalization.
    //!
    //! Specifically, the normalization ensures that:
    //! \f[ \max_{i=0}^n \left(\frac{|a_i|}{S}\right) = 1 \f]
    //! where  \f$ S \f$ is the maximum absolute value of
    //! the coefficients of the polynomial:
    //!
    //! \f[ S = \max_{i=0}^n |a_i| \f]
    //!
    //! If the maximum absolute value of the coefficients
    //! is zero (i.e., the polynomial is zero),
    //! the method will return zero without altering the coefficients.
    //!
    //! \return The scaling factor \f$ S \f$, which is used to normalize the
    //! polynomial coefficients.
    //!
    //! **Example**
    //!
    //! \code{cpp}
    //! Poly<double> P;
    //! P << 2, -3, 5; // P(x) = 2 - 3x + 5x^2
    //! Real scale = P.normalize(); // Scale will be 5, and coefficients will be
    //! normalized
    //! // New coefficients will be: [0.4, -0.6, 1.0] (after normalization)
    //! \endcode
    //!
    Real normalize()
    {
      // search max module coeff
      if ( m_order > 0 )
      {
        Real S{ this->cwiseAbs().maxCoeff() };
        if ( S > 0 ) this->to_eigen() /= S;
        adjust_degree();
        return S;
      }
      else
      {
        return 1;
      }
    }

    //!
    //! Purges (sets to zero) the coefficients of the polynomial
    //! \f$ p(x) = \sum_{i=0}^n a_i x^i \f$ for which the absolute
    //! value is less than or equal to a specified threshold \f$ \epsilon \f$.
    //!
    //! This method modifies the polynomial by removing insignificant
    //! coefficients, which can be useful for simplifying the polynomial
    //! and improving numerical stability.
    //! Specifically, it sets:
    //! \f[ a_i = 0 \quad \text{if} \quad |a_i| \leq \epsilon \f]
    //!
    //! This operation can be particularly beneficial when dealing
    //! with polynomials that have been generated through numerical methods,
    //! where small coefficients may result from rounding errors
    //! or other numerical inaccuracies.
    //!
    //! \param epsi The threshold value below which coefficients will be purged.
    //!
    //! **Example**
    //!
    //! \code{cpp}
    //! Poly<double> P;
    //! P << 0.1, 0.0005, -0.2, 0.3, 0.0001; // P(x) = 0.1 + 0.0005x - 0.2x^2 +
    //! 0.3x^3 + 0.0001x^4 P.purge(0.001); // Coefficients less than or equal to
    //! 0.001 will be set to 0
    //! // Resulting coefficients will be: [0.1, 0, -0.2, 0.3, 0]
    //! \endcode
    //!
    void purge( Real epsi )
    {
      if ( m_order > 0 )
      {
        Real MX = this->cwiseAbs().maxCoeff();
        if ( MX < 1 ) MX = 1;
        Real EPS = epsi * MX;
        for ( Integer i = 0; i < m_order; ++i )
        {
          Real & ai{ this->coeffRef( i ) };
          if ( std::abs( ai ) <= EPS ) ai = 0;
        }
      }
      adjust_degree();
    }

    //!
    //! Adjusts the polynomial order of \f$ p(x) = \sum_{i=0}^n a_i x^i \f$
    //! such that the leading coefficient \f$ a_n \f$ is non-zero.
    //! This method ensures that the polynomial's degree
    //! accurately reflects its meaningful terms by removing
    //! any trailing zero coefficients.
    //!
    //! Specifically, the method checks for any leading zeros in
    //! the coefficient array and modifies the polynomial's order
    //! to the highest degree where the coefficient is non-zero.
    //! This is crucial for maintaining the integrity of the polynomial
    //! representation and ensuring that operations involving the polynomial
    //! (such as evaluation or differentiation) yield correct results.
    //!
    //! If all coefficients are zero, the polynomial order will be adjusted to
    //! reflect this.
    //!
    //! **Example**
    //!
    //! \code{cpp}
    //! Poly<double> P;
    //! P.set_degree(4)
    //! P << 0, 0, 3, 4, 0; // P(x) = 3x^2 + 4x^3
    //! P.adjust_degree(); // Adjusts degree, ensuring a_n is non-zero
    //! // Resulting coefficients will reflect the highest degree with non-zero
    //! coefficient
    //! \endcode
    //!
    void adjust_degree()
    {
      while ( m_order > 0 && this->coeff( m_order - 1 ) == 0 ) --m_order;
      this->conservativeResize( m_order );
    }

    //!
    //! Counts the number of sign variations in the coefficients
    //! of the polynomial \f$ p(x) = \sum_{i=0}^n a_i x^i \f$.
    //! This method computes the number of sign changes in the
    //! sequence of coefficients \f$ [a_0, a_1, \ldots, a_n] \f$.
    //!
    //! This method returns the total count of sign changes in
    //! the coefficient sequence, which can be useful for
    //! root-finding algorithms and analyzing the behavior of the polynomial.
    //!
    //! \return The number of sign variations in the polynomial coefficients.
    //!
    //! **Example**
    //!
    //! \code{cpp}
    //! Poly<double> P;
    //! P << 3, -1, 2, -4; // P(x) = 3 - x + 2x^2 - 4x^3
    //! Integer variations = P.sign_variations(); // Counts the sign variations
    //! in [3, -1, 2, -4]
    //! // variations will be 3, corresponding to the sign changes: + -> - -> +
    //! -> -
    //! \endcode
    //!
    Integer sign_variations() const
    {
      Integer sign_var{ 0 };
      Integer last_sign{ 0 };
      for ( Integer i = 0; i < m_order; ++i )
      {
        Real v = this->coeff( i );
        if ( v > 0 )
        {
          if ( last_sign == -1 ) ++sign_var;
          last_sign = 1;
        }
        else if ( v < 0 )
        {
          if ( last_sign == 1 ) ++sign_var;
          last_sign = -1;
        }
      }
      return sign_var;
    }

    //!
    //! Transforms the polynomial \f$ p(x) = \sum_{i=0}^{n} a_i x^i \f$
    //! into a monic polynomial, which is defined as a polynomial
    //! where the leading coefficient is equal to 1.
    //! The method modifies the polynomial such that it takes the form:
    //! \f[ p(x) = x^n + \sum_{i=0}^{n-1} a_i x^i \f]
    //!
    //! Specifically, this method performs the following operations:
    //!
    //! 1. Divides all coefficients of the polynomial by the current
    //!    leading coefficient (i.e., \f$ a_n \f$).
    //!
    //! 2. Sets the leading coefficient (i.e., the coefficient of \f$ x^n \f$)
    //! to 1.
    //!
    //! **Example**
    //!
    //! \code{cpp}
    //! Poly<double> P;
    //! P << 4, -2, 1; // P(x) = 4 + (-2)x + 1x^2
    //! P.make_monic(); // Transforms P into monic form
    //! // Resulting coefficients will be: [1, -0.5, 0.25] corresponding to x^2
    //! - 0.5x + 0.25
    //! \endcode
    //!
    void make_monic()
    {
      this->to_eigen() /= this->coeff( m_order - 1 );
      this->coeffRef( m_order - 1 ) = 1;
    }

    //!
    //! Assignment operator for the polynomial class.
    //! This method copies the values from another polynomial
    //! \f$ p(x) \f$ into the current instance, ensuring that
    //! the internal state is updated to match the source polynomial.
    //!
    //! \param c The polynomial to be copied.
    //!
    //! \return A reference to the current polynomial instance after the
    //! assignment.
    //!
    //! **Example**
    //!
    //! \code{cpp}
    //! Poly<double> P, Q;
    //! P << 1, -2, 3; // P(x) = 1 - 2x + 3x^2
    //! Q = P; // Assigns P to Q
    //! // Now Q(x) = 1 - 2x + 3x^2
    //! \endcode
    //!
    Poly_t & operator=( Poly_t const & c )
    {
      this->resize( c.m_order );
      this->to_eigen().noalias() = c.to_eigen();
      m_order                    = c.m_order;
      return *this;
    }

    //!
    //! Unary negation operator for the polynomial.
    //! This method returns a new polynomial that represents
    //! the negation of the current polynomial \f$ p(x) \f$.
    //! The resulting polynomial is given by:
    //! \f[ -p(x) = -\sum_{i=0}^{n} a_i x^i \f]
    //!
    //! \return A new polynomial that is the negation of the current polynomial.
    //!
    //! **Example**
    //!
    //! \code{cpp}
    //! Poly<double> P;
    //! P << 3, -1, 2; // P(x) = 3 - x + 2x^2
    //! Poly<double> NegP = -P; // Negates P
    //! // NegP(x) = -3 + x - 2x^2
    //! \endcode
    //!
    Poly_t operator-() const { return Poly( -this->to_eigen() ); }

    //!
    //! Addition assignment operator for the polynomial.
    //! This method adds another polynomial \f$ q(x) \f$ to the
    //! current polynomial \f$ p(x) \f$.
    //! The operation modifies the current polynomial as follows:
    //! \f[ p(x) \leftarrow p(x) + q(x) \f]
    //!
    //! \param q The polynomial to be added.
    //! \return A reference to the current polynomial instance after the
    //! addition.
    //!
    //! **Example**
    //!
    //! \code{cpp}
    //! Poly<double> P, Q;
    //! P << 1, 2; // P(x) = 1 + 2x
    //! Q << 3, -4; // Q(x) = 3 - 4x
    //! P += Q; // Adds Q to P
    //! // Now P(x) = 4 - 2x
    //! \endcode
    //!
    Poly_t & operator+=( Poly_t const & q )
    {
      Integer max_order = std::max( m_order, q.m_order );
      Integer min_order = std::min( m_order, q.m_order );

      // ridimensiona vettore coefficienti senza distruggere il contenuto
      this->conservativeResize( max_order );

      // somma i coefficienti fino al grado comune ad entrambi i polinomi
      this->head( min_order ).noalias() += q.head( min_order );

      if ( Integer n_tail{ q.m_order - m_order }; n_tail > 0 ) this->tail( n_tail ).noalias() = q.tail( n_tail );

      m_order = max_order;
      return *this;
    }

    //!
    //! Subtraction assignment operator for the polynomial.
    //! This method subtracts another polynomial \f$ q(x) \f$ from the current
    //! polynomial \f$ p(x) \f$. The operation modifies the current polynomial
    //! as follows:
    //! \f[ p(x) \leftarrow p(x) - q(x) \f]
    //!
    //! \param q The polynomial to be subtracted.
    //!
    //! \return A reference to the current polynomial instance after the
    //! subtraction.
    //!
    //! **Example**
    //!
    //! \code{cpp}
    //! Poly<double> P, Q;
    //! P << 5, 0, 2; // P(x) = 5 + 0*x + 2*x^2
    //! Q << 1, 1; // Q(x) = 1 + x
    //! P -= Q; // Subtracts Q from P
    //! // Now P(x) = 4 - x + 2*x^2
    //! \endcode
    //!
    Poly_t & operator-=( Poly_t const & q )
    {
      Integer max_order{ std::max( m_order, q.m_order ) };
      Integer min_order{ std::min( m_order, q.m_order ) };

      // ridimensiona vettore coefficienti senza distruggere il contenuto
      this->conservativeResize( max_order );

      // somma i coefficienti fino al grado comune ad entrambi i polinomi
      this->head( min_order ).noalias() -= q.head( min_order );

      if ( Integer n_tail{ q.m_order - m_order }; n_tail > 0 ) this->tail( n_tail ).noalias() = -q.tail( n_tail );

      m_order = max_order;
      return *this;
    }

    //!
    //! Multiplication assignment operator for the polynomial.
    //! This method multiplies the current polynomial \f$ p(x) \f$
    //! by another polynomial \f$ q(x) \f$.
    //! The operation modifies the current polynomial as follows:
    //! \f[ p(x) \leftarrow p(x) \cdot q(x) \f]
    //!
    //! \param q The polynomial to be multiplied.
    //!
    //! \return A reference to the current polynomial instance after the
    //! multiplication.
    //!
    //! **Example**
    //!
    //! \code{cpp}
    //! Poly<double> P, Q;
    //! P << 2, 1; // P(x) = 2 + x
    //! Q << 1, -1; // Q(x) = 1 - x
    //! P *= Q; // Multiplies P by Q
    //! // Now P(x) = (2 + x)(1 - x) = 2 - 2*x + x^2
    //! \endcode
    //!
    Poly_t & operator*=( Poly_t const & q )
    {
      dvec_t        a( this->to_eigen() );  // fa una copia dei coefficienti del vettore
      Integer const new_order{ m_order + q.m_order - 1 };
      this->resize( m_order + q.m_order - 1 );  // nuovo polinomio contenente il risultato
      this->setZero();
      for ( Integer i = 0; i < m_order; ++i )
        for ( Integer j = 0; j < q.m_order; ++j ) this->coeffRef( i + j ) += a.coeff( i ) * q.coeff( j );
      m_order = new_order;
      return *this;
    }

    //!
    //! Addition assignment operator for the polynomial with a scalar.
    //! This method adds a scalar value \f$ a \f$ to the current polynomial \f$
    //! p(x) \f$. The operation modifies the current polynomial as follows:
    //! \f[ p(x) \leftarrow p(x) + a \f]
    //!
    //! \param a The scalar value to be added.
    //!
    //! \return A reference to the current polynomial instance after the
    //! addition.
    //!
    //! **Example**
    //!
    //! \code{cpp}
    //! Poly<double> P;
    //! P << 1, 2; // P(x) = 1 + 2x
    //! P += 3; // Adds 3 to P
    //! // Now P(x) = 4 + 2x
    //! \endcode
    //!
    Poly_t & operator+=( Real a )
    {
      if ( m_order > 0 )
        this->coeffRef( 0 ) += a;
      else
      {
        this->resize( 1 );
        this->coeffRef( 0 ) = a;
        m_order             = 1;
      }
      return *this;
    }

    //!
    //! Subtraction assignment operator for the polynomial with a scalar.
    //! This method subtracts a scalar value \f$ a \f$ from the
    //! current polynomial \f$ p(x) \f$.
    //! The operation modifies the current polynomial as follows:
    //! \f[ p(x) \leftarrow p(x) - a \f]
    //!
    //! \param a The scalar value to be subtracted.
    //!
    //! \return A reference to the current polynomial instance after the
    //! subtraction.
    //!
    //! **Example**
    //!
    //! \code{cpp}
    //! Poly<double> P;
    //! P << 4, -1; // P(x) = 4 - x
    //! P -= 2; // Subtracts 2 from P
    //! // Now P(x) = 2 - x
    //! \endcode
    //!
    Poly_t & operator-=( Real a )
    {
      if ( m_order > 0 )
        this->coeffRef( 0 ) -= a;
      else
      {
        this->resize( 1 );
        this->coeffRef( 0 ) = -a;
        m_order             = 1;
      }
      return *this;
    }

    //!
    //! Multiplication assignment operator for the polynomial with a scalar.
    //! This method multiplies the current polynomial \f$ p(x) \f$
    //! by a scalar value \f$ a \f$.
    //! The operation modifies the current polynomial as follows:
    //! \f[ p(x) \leftarrow p(x) \cdot a \f]
    //!
    //! \param a The scalar value to be multiplied.
    //!
    //! \return A reference to the current polynomial instance after the
    //! multiplication.
    //!
    //! **Example**
    //!
    //! \code{cpp}
    //! Poly<double> P;
    //! P << 1, 2; // P(x) = 1 + 2x
    //! P *= 2; // Multiplies P by 2
    //! // Now P(x) = 2 + 4x
    //! \endcode
    //!
    Poly_t & operator*=( Real a )
    {
      this->to_eigen() *= a;
      return *this;
    }
  };

  //!
  //! \ingroup Zeros
  //!
  //! \brief Class for managing and computing the Sturm sequence associated with
  //! a polynomial.
  //!
  //! The Sturm class is dedicated to the construction and manipulation of the
  //! Sturm sequence for a given polynomial \f$ p(x) \f$. This sequence is
  //! instrumental in real algebraic geometry, as it provides crucial insights
  //! into the behavior of polynomial roots. Specifically, it enables users to
  //! ascertain the number of real roots within specific intervals and helps
  //! identify the locations of distinct real roots [2,3].
  //!
  //! \tparam Real The data type used for the polynomial coefficients, typically
  //! a floating-point type.
  //!
  //! **Features**
  //!
  //! - **Sturm Sequence Construction**: Efficiently builds the Sturm sequence
  //! for a specified polynomial, facilitating root analysis.
  //!
  //! - **Sign Variation Counting**: Provides methods to compute the number of
  //! sign variations of the Sturm sequence at given points,
  //!   which is essential for determining the count of real roots in specified
  //!   intervals.
  //!
  //! - **Root Isolation**: Implements techniques to isolate intervals that
  //! contain single roots, aiding in numerical root-finding methods [4].
  //!
  //! - **Interval Analysis**: Allows for the examination of polynomial behavior
  //! across intervals, enhancing the understanding of root distribution.
  //!
  template <typename Real> class Sturm
  {
  public:
    using Integer = int;
    using Poly_t  = Poly<Real>;
    using dvec_t  = Eigen::Matrix<Real, Eigen::Dynamic, 1>;

    //! Structure representing an interval with sign variation counts
    using Interval = struct Interval
    {
      Real    a;          ///< Left endpoint
      Real    b;          ///< Right endpoint
      Integer va;         ///< Sign variations at left endpoint
      Integer vb;         ///< Sign variations at right endpoint
      bool    a_on_root;  ///< True if left endpoint is a root
      bool    b_on_root;  ///< True if right endpoint is a root
    };

    //! Function wrapper for bracket-based root refinement
    class Bracket_fun : public Bracket_base_fun<Real>
    {
      Poly<Real> const * P = nullptr;

    public:
      void setup( Poly<Real> const * Pin ) { P = Pin; }
      Real eval( Real x ) const override { return P->eval( x ); }
    };

  private:
    AlgoBracket<Real> m_solver;
    Bracket_fun       m_fun;

    vector<Poly_t>   m_sturm;      ///< Sturm sequence polynomials
    vector<Interval> m_intervals;  ///< Isolated intervals containing roots
    dvec_t           m_roots;      ///< Refined roots
    Real             m_a{ 0 };     ///< Left bound of search interval
    Real             m_b{ 0 };     ///< Right bound of search interval

  public:
    Sturm() = default;

    //!
    //! \brief Constructs the Sturm sequence for a given polynomial.
    //!
    //! This method takes a polynomial \f$ p(x) \f$ as input and builds its
    //! corresponding Sturm sequence using the Euclidean algorithm [2,3].
    //!
    //! \param p The polynomial \f$ p(x) \f$ for which to build the Sturm
    //! sequence.
    //!          This should be a valid polynomial of type Poly_t.
    //!
    //! \note The Sturm sequence will be stored internally within the class,
    //!       allowing subsequent methods to utilize it for root counting
    //!       or determining the nature of the roots of the polynomial.
    //!
    //! \warning Ensure that the polynomial is non-constant and valid,
    //!          as the Sturm sequence is not defined for constant polynomials.
    //!
    void build( Poly_t const & p )
    {
      m_intervals.clear();
      Poly_t DP, M, R;
      p.derivative( DP );
      m_sturm.clear();
      m_sturm.reserve( p.order() );
      m_sturm.emplace_back( p );
      m_sturm.back().adjust_degree();
      m_sturm.emplace_back( DP );
      m_sturm.back().adjust_degree();
      Integer ns{ 1 };
      while ( true )
      {
        divide( m_sturm[ns - 1], m_sturm[ns], M, R );
        if ( R.order() <= 0 ) break;
        m_sturm.emplace_back( -R );
        ++ns;
      }
      // divide by GCD
      for ( Integer i{ 0 }; i < ns; ++i )
      {
        divide( m_sturm[i], m_sturm.back(), M, R );
        M.normalize();
        m_sturm[i] = M;
      }
      m_sturm.back().set_scalar( 1 );
    }

    //!
    //! \brief Retrieves the length of the stored Sturm sequence.
    //!
    //! This method returns the number of polynomials in the Sturm sequence
    //! that has been built for the corresponding polynomial.
    //!
    //! \return The number of polynomials in the Sturm sequence as an
    //!         integer value.
    //!
    //! \note This method does not modify the Sturm sequence or its
    //!       contents. It only provides information about its size.
    //!
    Integer length() const { return Integer( m_sturm.size() ); }

    //!
    //! \brief Retrieves the i-th polynomial of the stored Sturm sequence.
    //!
    //! This method returns a constant reference to the polynomial at the
    //! specified index `i` in the Sturm sequence. The index is zero-based,
    //! meaning that `get(0)` will return the first polynomial in the
    //! sequence, `get(1)` will return the second, and so on.
    //!
    //! \param i The index of the polynomial to retrieve from the Sturm
    //! sequence.
    //!          It must be a valid index (0 <= i < length of Sturm sequence).
    //!
    //! \return A constant reference to the i-th polynomial in the Sturm
    //! sequence.
    //!
    //! \throws std::out_of_range if the index `i` is out of bounds.
    //!
    //! \note This method does not modify the Sturm sequence or its contents.
    //!
    //! **Example**
    //!
    //! \code
    //! Sturm<Poly_t> sturm;
    //! Poly_t polynomial = ...; // Assume polynomial is initialized
    //! sturm.build(polynomial); // Build the Sturm sequence from the polynomial
    //!
    //! // Get the first polynomial in the Sturm sequence
    //! Poly_t const & first_poly = sturm.get(0);
    //!
    //! // Print or utilize the first polynomial
    //! std::cout << "First polynomial: " << first_poly << '\n';
    //! \endcode
    //!
    Poly_t const & get( Integer i ) const { return m_sturm[i]; }

    //!
    //! \brief Computes the sign variations of the stored Sturm sequence at a
    //! given point \f$ x \f$.
    //!
    //! This method evaluates the sign variations of the polynomials in the
    //! stored Sturm sequence at the specified point \f$ x \f$. The number of
    //! sign variations can help determine the number of real roots of the
    //! polynomial in the interval up to \f$ x \f$ (Sturm's theorem [2]).
    //!
    //! \param x The point at which to compute the sign variations of the Sturm
    //! sequence.
    //! \param on_root A reference to a boolean variable that will be set to
    //! true if \f$ x \f$
    //!                is a root of the first polynomial of the Sturm sequence;
    //!                otherwise, it will be false.
    //!
    //! \return The number of sign variations in the Sturm sequence at the point
    //! \f$ x \f$.
    //!
    //! **Example**
    //!
    //! \code
    //! Sturm<Poly_t> sturm;
    //! Poly_t polynomial = ...; // Assume polynomial is initialized
    //! sturm.build(polynomial); // Build the Sturm sequence from the polynomial
    //!
    //! Real x = ...; // The point at which to evaluate the sign variations
    //! bool is_on_root = false;
    //!
    //! Integer variations = sturm.sign_variations(x, is_on_root);
    //!
    //! std::cout << "Sign variations at x = " << x << ": " << variations <<
    //! '\n'; if (is_on_root) {
    //!   std::cout << "The point x is a root of polynomial.\n";
    //! } else {
    //!   std::cout << "The point x is not a root of polynomial.\n";
    //! }
    //! \endcode
    Integer sign_variations( Real x, bool & on_root ) const
    {
      Integer const npoly{ static_cast<Integer>( m_sturm.size() ) };
      Integer       sign_var{ 0 };
      Integer       last_sign{ 0 };
      Real          v{ m_sturm[0].eval( x ) };
      on_root = false;
      if ( v > 0 )
        last_sign = 1;
      else if ( v < 0 )
        last_sign = -1;
      else
      {
        on_root   = true;
        last_sign = 0;
      }
      for ( Integer i{ 1 }; i < npoly; ++i )
      {
        v = m_sturm[i].eval( x );
        if ( v > 0 )
        {
          if ( last_sign == -1 ) ++sign_var;
          last_sign = 1;
        }
        else if ( v < 0 )
        {
          if ( last_sign == 1 ) ++sign_var;
          last_sign = -1;
        }
      }
      return sign_var;
    }

    //!
    //! \brief Computes subintervals containing single roots within a given
    //! interval \f$ [a,b] \f$.
    //!
    //! This method identifies subintervals within the specified range \f$ [a,b]
    //! \f$ that each contain a single root of the polynomial. It performs root
    //! isolation using Sturm's theorem and bisection, and returns the number of
    //! such intervals found.
    //!
    //! \param a The lower bound of the interval.
    //! \param b The upper bound of the interval.
    //!
    //! \return The number of intervals (roots) found within the interval \f$
    //! [a,b] \f$.
    //!
    //! **example**
    //!
    //! \code
    //! Sturm<Poly_t> sturm;
    //! Poly_t polynomial = ...; // Assume polynomial is initialized
    //! sturm.build(polynomial); // Build the Sturm sequence from the polynomial
    //!
    //! Real a = -1.0; // Lower bound of the interval
    //! Real b = 1.0;  // Upper bound of the interval
    //!
    //! Integer root_count = sturm.separate_roots(a, b);
    //! std::cout << "Number of roots found in the interval [" << a << ", " << b
    //! << "]: " << root_count << '\n';
    //! \endcode
    //!
    Integer separate_roots( Real a, Real b )
    {
      using std::abs;
      using std::max;

      m_intervals.clear();
      m_intervals.reserve( m_sturm.size() );

      Interval I0, I1;
      m_a = I0.a = a;
      m_b = I0.b = b;

      I0.va = sign_variations( I0.a, I0.a_on_root );
      I0.vb = sign_variations( I0.b, I0.b_on_root );

      Integer n_roots{ std::abs( I0.va - I0.vb ) };

      if ( n_roots <= 1 )
      {
        if ( n_roots == 1 && !I0.a_on_root && !I0.b_on_root ) { m_intervals.push_back( I0 ); }
        if ( I0.a_on_root )
        {
          I1.a = I1.b = I0.a;
          I1.va = I1.vb = I0.va;
          I1.a_on_root = I1.b_on_root = true;
          m_intervals.push_back( I1 );
        }
        if ( I0.b_on_root )
        {
          I1.a = I1.b = I0.b;
          I1.va = I1.vb = I0.vb;
          I1.a_on_root = I1.b_on_root = true;
          m_intervals.push_back( I1 );
        }
        return static_cast<Integer>( m_intervals.size() );
      }

      // search intervals
      vector<Interval> I_stack;
      I_stack.clear();
      I_stack.reserve( m_sturm.size() );
      I_stack.push_back( I0 );
      while ( I_stack.size() > 0 )
      {
        I0 = I_stack.back();
        I_stack.pop_back();
        // controllo se una sola radice
        n_roots = std::abs( I0.va - I0.vb );
        if ( n_roots <= 1 )
        {
          if ( I0.a_on_root )
          {
            I0.b         = I0.a;
            I0.vb        = I0.va;
            I0.b_on_root = true;
            m_intervals.push_back( I0 );
          }
          else if ( I0.b_on_root )
          {
            I0.a         = I0.b;
            I0.va        = I0.vb;
            I0.a_on_root = true;
            m_intervals.push_back( I0 );
          }
          else if ( n_roots == 1 ) { m_intervals.push_back( I0 ); }
        }
        else if ( abs( I0.b - I0.a ) <= 10 * machine_eps<Real>() * max( Real( 1 ), max( abs( I0.b ), abs( I0.a ) ) ) )
        {
          I1.a = I1.b = I0.a;
          I1.va = I1.vb = 0;
          I1.a_on_root = I1.b_on_root = true;
          I_stack.push_back( I1 );
        }
        else
        {
          Real    c{ ( I0.a + I0.b ) / 2 };
          bool    c_on_root;
          Integer vc = sign_variations( c, c_on_root );
          // check interval [a,c]
          if ( I0.va != vc || c_on_root || I0.a_on_root )
          {
            if ( c < I0.b )
            {  // check if it is a true reduction
              I1.a         = I0.a;
              I1.va        = I0.va;
              I1.a_on_root = I0.a_on_root;
              I1.b         = c;
              I1.vb        = vc;
              I1.b_on_root = c_on_root;
              I_stack.push_back( I1 );
            }
            else if ( c_on_root )
            {
              I1.a = I1.b  = c;
              I1.a_on_root = I1.b_on_root = true;
              I1.va = I1.vb = 0;
              I_stack.push_back( I1 );
            }
            else if ( I0.a_on_root )
            {
              I1.a = I1.b  = I0.a;
              I1.a_on_root = I1.b_on_root = true;
              I1.va = I1.vb = 0;
              I_stack.push_back( I1 );
            }
          }
          // check interval [c,b]
          if ( I0.vb != vc || I0.b_on_root )
          {
            if ( c > I0.a )
            {
              I1.a         = c;
              I1.va        = vc;
              I1.a_on_root = c_on_root;
              I1.b         = I0.b;
              I1.vb        = I0.vb;
              I1.b_on_root = I0.b_on_root;
              I_stack.push_back( I1 );
            }
            else if ( I0.b_on_root )
            {
              I1.a = I1.b  = I0.b;
              I1.a_on_root = I1.b_on_root = true;
              I1.va = I1.vb = 0;
              I_stack.push_back( I1 );
            }
          }
        }
      }
      // sort intervals
      std::sort(
        m_intervals.begin(),
        m_intervals.end(),
        []( Interval const & Sa, Interval const & Sb ) { return Sa.a < Sb.a; } );
      return static_cast<Integer>( m_intervals.size() );
    }

    //!
    //! \brief Compute an interval \f$ [a,b] \f$ that contains all real roots
    //! and separate them into subintervals.
    //!
    //! This method uses Cauchy's bounds to compute an interval \f$ [-B, B] \f$,
    //! where \f$ B \f$ is an upper bound for the absolute values of the real
    //! roots of the polynomial. It then isolates each root by dividing this
    //! interval into subintervals, each containing exactly one root. It returns
    //! the number of subintervals (i.e., roots) found.
    //!
    //! The method leverages Cauchy's bound \f$ B = 1 + \frac{\max |a_i|}{|a_n|}
    //! \f$, where \f$ a_n \f$ is the leading coefficient and \f$ a_i \f$ are
    //! the other coefficients, to determine an interval that contains all real
    //! roots.
    //!
    //! \return The number of subintervals (roots) found.
    //!
    //! **Example**
    //!
    //! \code
    //! Sturm<Poly_t> sturm;
    //! Poly_t polynomial = ...; // Assume polynomial is initialized
    //! sturm.build(polynomial); // Build the Sturm sequence from the polynomial
    //!
    //! Integer root_count = sturm.separate_roots(); // Automatically separate
    //! all real roots within computed bounds std::cout << "Number of real roots
    //! found: " << root_count << '\n';
    //! \endcode
    //!
    //! This method can be used to isolate all real roots of a polynomial in a
    //! simple, automated fashion.
    //!
    Integer separate_roots()
    {
      // Cauchy's bounds for roots
      Real an  = m_sturm[0].leading_coeff();
      Real bnd = 1 + m_sturm[0].cwiseAbs().maxCoeff() / abs( an );
      return separate_roots( -bnd, bnd );
    }

    //!
    //! \brief Returns the number of roots found by the Sturm sequence.
    //!
    //! After running root separation (e.g., using `separate_roots()`), this
    //! method returns the number of intervals found, each containing a single
    //! root.
    //!
    //! \return The number of root-containing intervals.
    //!
    //! **Example**
    //!
    //! \code
    //! Sturm<Poly_t> sturm;
    //! Poly_t polynomial = ...; // Initialize polynomial
    //! sturm.build(polynomial);
    //! sturm.separate_roots();
    //!
    //! Integer num_roots = sturm.n_roots();
    //! std::cout << "Number of real roots: " << num_roots << '\n';
    //! \endcode
    //!
    Integer n_roots() const { return Integer( m_intervals.size() ); }

    //!
    //! \brief Returns the left boundary \f$ a \f$ of the interval \f$ [a,b] \f$
    //! containing all real roots.
    //!
    //! This method returns the left endpoint of the interval that contains all
    //! real roots of the polynomial after using methods such as
    //! `separate_roots()`.
    //!
    //! \return The left boundary \f$ a \f$ of the interval.
    //!
    Real a() const { return m_a; }

    //!
    //! \brief Returns the right boundary \f$ b \f$ of the interval \f$ [a,b]
    //! \f$ containing all real roots.
    //!
    //! This method returns the right endpoint of the interval that contains all
    //! real roots of the polynomial after using methods such as
    //! `separate_roots()`.
    //!
    //! \return The right boundary \f$ b \f$ of the interval.
    //!
    Real b() const { return m_b; }

    //!
    //! \brief Returns the i-th interval containing a single root.
    //!
    //! After separating the roots of the polynomial, this method provides
    //! access to the specific interval that contains the i-th root. The
    //! interval is defined by a structure `Interval`, which holds the
    //! boundaries of the interval.
    //!
    //! \param i The index of the interval (root).
    //! \return The interval \f$ [a_i, b_i] \f$ containing the i-th root.
    //!
    Interval const & get_interval( Integer i ) const { return m_intervals[i]; }

    //!
    //! \brief Refine the roots of the polynomial after the intervals have been
    //! separated.
    //!
    //! This method computes the roots within the intervals determined by the
    //! `separate_roots()` method. After separating the intervals that contain
    //! the roots, this method applies a refinement procedure to accurately
    //! locate the roots within those intervals using a bracket-based root
    //! finding algorithm.
    //!
    //! **Example**
    //!
    //! \code
    //! Sturm<Poly_t> sturm;
    //! Poly_t polynomial = ...; // Initialize polynomial
    //! sturm.build(polynomial);
    //! sturm.separate_roots();
    //! sturm.refine_roots();
    //!
    //! // Now the roots are refined and can be accessed
    //! dvec_t const& refined_roots = sturm.roots();
    //! \endcode
    //!
    void refine_roots()
    {
      m_fun.setup( &m_sturm[0] );
      m_roots.resize( m_intervals.size() );
      Integer n{ 0 };
      for ( auto & I : m_intervals )
      {
        Real & r{ m_roots.coeffRef( n++ ) };
        if ( I.a_on_root )
          r = I.a;
        else if ( I.b_on_root )
          r = I.b;
        else
        {
          r = m_solver.eval( I.a, I.b, &m_fun );
          if ( !m_solver.converged() ) fmt::print( "Warning: Sturm<Real>::refine_roots failed at interval N.{}\n", n );
        }
      }
    }

    //!
    //! \brief Returns a vector containing the computed roots after refinement.
    //!
    //! This method returns the refined roots of the polynomial that were
    //! computed using the `refine_roots()` method. The roots are stored in a
    //! vector of type `dvec_t`.
    //!
    //! \return A constant reference to the vector containing the refined roots.
    //!
    //! **Example**
    //!
    //! \code
    //! // After separating and refining roots
    //! dvec_t const& roots = sturm.roots();
    //! for (const auto& root : roots) {
    //!     std::cout << "Root: " << root << '\n';
    //! }
    //! \endcode
    //!
    dvec_t const & roots() const { return m_roots; }
  };

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real> Poly<Real> operator+( Poly<Real> const & a, Poly<Real> const & b )
  {
    using Integer = typename Poly<Real>::Integer;
    Integer    max_order{ std::max( a.order(), b.order() ) };
    Integer    min_order{ std::min( a.order(), b.order() ) };
    Poly<Real> sum( max_order );  // nuovo polinomio contenente la somma

    // somma i coefficienti fino al grado comune ad entrambi i polinomi
    sum.head( min_order ).noalias() = a.head( min_order ) + b.head( min_order );
    Integer n_tail{ max_order - min_order };
    if ( n_tail > 0 )
    {
      if ( a.order() > b.order() )
        sum.tail( n_tail ).noalias() = a.tail( n_tail );
      else
        sum.tail( n_tail ).noalias() = b.tail( n_tail );
    }
    return sum;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real> Poly<Real> operator+( Poly<Real> const & a, Real b )
  {
    using Integer = typename Poly<Real>::Integer;
    Integer    max_order{ std::max( a.order(), 1 ) };
    Poly<Real> sum( max_order );  // nuovo polinomio contenente la somma

    // somma i coefficienti fino al grado comune ad entrambi i polinomi
    if ( a.order() > 0 )
    {
      sum.coeffRef( 0 ) = a.coeff( 0 ) + b;
      if ( a.order() > 1 ) sum.tail( a.order() - 1 ).noalias() = a.tail( a.order() - 1 );
    }
    else
    {
      sum.coeffRef( 0 ) = b;
    }
    return sum;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real> Poly<Real> operator+( Real a, Poly<Real> const & b )
  {
    using Integer = typename Poly<Real>::Integer;
    Integer    max_order{ std::max( b.order(), 1 ) };
    Poly<Real> sum( max_order );  // nuovo polinomio contenente la somma

    // somma i coefficienti fino al grado comune ad entrambi i polinomi
    if ( b.order() > 0 )
    {
      sum.coeffRef( 0 ) = a + b.coeff( 0 );
      if ( b.order() > 1 ) sum.tail( b.order() - 1 ).noalias() = b.tail( b.order() - 1 );
    }
    else
    {
      sum.coeffRef( 0 ) = a;
    }
    return sum;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real> Poly<Real> operator-( Poly<Real> const & a, Poly<Real> const & b )
  {
    using Integer = typename Poly<Real>::Integer;
    Integer    max_order{ std::max( a.order(), b.order() ) };
    Integer    min_order{ std::min( a.order(), b.order() ) };
    Poly<Real> sum( max_order );  // nuovo polinomio contenente la somma

    // somma i coefficienti fino al grado comune ad entrambi i polinomi
    sum.head( min_order ).noalias() = a.head( min_order ) - b.head( min_order );
    Integer n_tail                  = max_order - min_order;
    if ( n_tail > 0 )
    {
      if ( a.order() > b.order() )
        sum.tail( n_tail ).noalias() = a.tail( n_tail );
      else
        sum.tail( n_tail ).noalias() = -b.tail( n_tail );
    }
    return sum;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real> Poly<Real> operator-( Poly<Real> const & a, Real b )
  {
    using Integer = typename Poly<Real>::Integer;
    Integer    max_order{ std::max( a.order(), 1 ) };
    Poly<Real> sum( max_order );  // nuovo polinomio contenente la somma

    // somma i coefficienti fino al grado comune ad entrambi i polinomi
    if ( a.order() > 0 )
    {
      sum.coeffRef( 0 ) = a.coeff( 0 ) - b;
      if ( a.order() > 1 ) sum.tail( a.order() - 1 ).noalias() = a.tail( a.order() - 1 );
    }
    else
    {
      sum.coeffRef( 0 ) = -b;
    }
    return sum;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real> Poly<Real> operator-( Real a, Poly<Real> const & b )
  {
    using Integer = typename Poly<Real>::Integer;
    Integer    max_order{ std::max( b.order(), 1 ) };
    Poly<Real> sum( max_order );  // nuovo polinomio contenente la somma

    // somma i coefficienti fino al grado comune ad entrambi i polinomi
    if ( b.order() > 0 )
    {
      sum.coeffRef( 0 ) = a - b.coeff( 0 );
      if ( b.order() > 1 ) sum.tail( b.order() - 1 ).noalias() = -b.tail( b.order() - 1 );
    }
    else
    {
      sum.coeffRef( 0 ) = a;
    }
    return sum;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real> Poly<Real> operator*( Poly<Real> const & a, Poly<Real> const & b )
  {
    using Integer = typename Poly<Real>::Integer;
    Poly<Real> prd( a.order() + b.order() - 1 );  // nuovo polinomio contenente il risultato
    for ( Integer i{ 0 }; i < a.order(); ++i )
      for ( Integer j{ 0 }; j < b.order(); ++j ) prd.coeffRef( i + j ) += a.coeff( i ) * b.coeff( j );
    return prd;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real> Poly<Real> operator*( Real a, Poly<Real> const & b )
  {
    Poly<Real> prd( b.order() );  // nuovo polinomio contenente il risultato
    prd.noalias() = a * b.to_eigen();
    return prd;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real> Poly<Real> operator*( Poly<Real> const & a, Real b )
  {
    Poly<Real> prd( a.order() );  // nuovo polinomio contenente il risultato
    prd.noalias() = a.to_eigen() * b;
    return prd;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */

  //!
  //! Divide the polynomial \f$ p(x) \f$ by \f$ q(x) \f$
  //! with remainder, i.e. \f$ p(x) = s(x)q(x)+r(x) \f$.
  //!
  //! @tparam Real Floating-point type
  //! @param p Dividend polynomial
  //! @param q Divisor polynomial
  //! @param M Output: quotient polynomial
  //! @param R Output: remainder polynomial
  //!
  template <typename Real> void divide( Poly<Real> const & p, Poly<Real> const & q, Poly<Real> & M, Poly<Real> & R )
  {
    UTILS_ASSERT(
      p.order() > 0 && q.order(),
      "Poly::divide( p, q, N, R ) ∂p = {}, ∂q = {} must be greather of 0",
      p.order(),
      q.order() );

    using Integer = typename Poly<Real>::Integer;

    Poly<Real> P( p ), Q( q );
    //
    //  scale polynomials
    //
    //  P(x) = p(x) / scaleP, Q(x) = q(x) / scaleQ
    //
    Real scaleP = P.normalize();
    Real scaleQ = Q.normalize();

    //
    // P(x) = Q(x) * M(x) + R(x)
    //
    R = P;
    Real    lcQ{ Q.leading_coeff() };
    Integer dd{ R.order() - Q.order() };
    if ( dd < 0 )
    {
      // P = Q +R
      M.set_scalar( 1 );
      R = P - Q;
    }
    else
    {
      Integer R_degree = R.degree();
      M.set_order( dd + 1 );

      assert( lcQ != 0 && "Poly::divide(p,q,M,R), leading coefficient of q(x) is 0!" );

      while ( dd >= 0 && R_degree >= 0 )
      {
        Real lcR{ R( R_degree ) };
        Real bf{ lcR / lcQ };
        M.coeffRef( dd ) = bf;
        R.segment( dd, Q.degree() ).noalias() -= bf * Q.head( Q.degree() );
        R.coeffRef( R_degree ) = 0;
        --R_degree;
        --dd;
      }

      // don not purge remainder
      // this can be done externally
      // R.purge(epsi);
      R.adjust_degree();
    }

    // scale back polinomials
    //
    // P(x) = Q(x) * M(x) + R(x)
    // p(x) / scaleP = q(x) / scaleQ * M(x) + R(x)
    // p(x) = q(x) * (scaleP/scaleQ) * M(x) + scaleP*R(x)
    //
    M *= scaleP / scaleQ;
    R *= scaleP;
  }

  //!
  //! Given \f$ p(x) \f$ and \f$ q(x) \f$ compute G.C.D, i.e.
  //! the polynomial \f$ g(x) \f$ such that \f$ q(x) | p(x) \f$
  //! and \f$ g(x) | q(x) \f$ and if another polynomial
  //! \f$ h(x) \f$ is such that \f$ h(x) | p(x) \f$
  //! and \f$ h(x) | q(x) \f$ then \f$ h(x) | g(x) \f$.
  //!
  //! @tparam Real Floating-point type
  //! @param p First polynomial
  //! @param q Second polynomial
  //! @param g Output: GCD polynomial
  //! @param epsi Tolerance for coefficient purging (default: 1e-20)
  //!
  template <typename Real>
  void GCD( Poly<Real> const & p, Poly<Real> const & q, Poly<Real> & g, Real epsi = Real( 1e-20 ) )
  {
    if ( q.order() > 0 )
    {
      Poly<Real> M, R;
      divide( p, q, M, R );
      R.purge( epsi );
      GCD( q, R, g, epsi );
    }
    else
    {
      g = p;
    }
    g.normalize();
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real> inline string Poly<Real>::to_string() const
  {
    if ( this->order() <= 0 ) return "EMPTY!";
    if ( this->order() == 1 ) return fmt::format( "{}", this->coeff( 0 ) );
    if ( ( *this ).cwiseAbs().maxCoeff() == 0 ) return "0";

    bool   empty = true;  // true indica che i coefficienti finora sono nulli
    string s     = "";    // segno
    Real   c     = 0;     // coefficiente
    string e     = "";    // esponente
    string res   = "";

    // controlla se esiste il primo coefficiente (grado 0)
    if ( this->coeff( 0 ) != 0 )
    {
      res   = fmt::format( "{}", this->coeff( 0 ) );
      empty = false;
    }

    for ( typename Poly<Real>::Integer i = 1; i < this->order(); ++i )
    {
      // se il coefficiente e` negativo...
      if ( this->coeff( i ) < 0 )
      {
        // e se i coefficienti precenti erano nulli...
        if ( empty )
        {
          s     = "";  // ...non scrive il segno
          c     = this->coeff( i );
          empty = false;
        }
        else
        {
          s = " - ";              // ...altrimenti scrive il segno come separatore
          c = -this->coeff( i );  // e inverte il segno del coefficiente
        }

        // se il coefficiente e` positivo...
      }
      else if ( this->coeff( i ) > 0 )
      {
        c = this->coeff( i );
        // e se i coefficienti precenti erano nulli...
        if ( empty )
        {
          s     = "";  // ...non scrive il segno
          empty = false;
        }
        else
        {
          s = " + ";  // ...altrimenti scrive il segno come separatore
        }

        // se il coefficiente e` zero...
      }
      else
      {
        continue;  // ...procede al prossimo
      }

      // se il grado e` 1 non scrive l'esponente
      if ( i == 1 )
        e = "x";
      else
        e = fmt::format( "x^{}", i );

      // se il coeff è 1 non lo stampo
      if ( c == 1 )
      {
        res += s;
        res += e;
      }
      else
        res += fmt::format( "{}{} {}", s, c, e );
    }
    return res;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real, typename Char>
  inline std::basic_ostream<Char> & operator<<( std::basic_ostream<Char> & output, Poly<Real> const & p )
  {
    output << p.to_string();
    return output;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real, typename Char>
  inline std::basic_ostream<Char> & operator<<( std::basic_ostream<Char> & output, Sturm<Real> const & S )
  {
    using Integer = typename Poly<Real>::Integer;
    output << "Sturm sequence\n";
    for ( Integer i = 0; i < S.length(); ++i ) output << "P_" << i << "(x) = " << S.get( i ) << '\n';

    Integer n = S.n_roots();
    if ( n > 0 )
    {
      fmt::print( output, "roots separation for interval [{},{}]\n", S.a(), S.b() );
      for ( Integer i = 0; i < n; ++i )
      {
        typename Sturm<Real>::Interval const & I = S.get_interval( i );
        fmt::print( output, "I = [{}, {}], V = [{}, {}]\n", I.a, I.b, I.va, I.vb );
      }
    }
    return output;
  }

}  // namespace Utils

namespace fmt
{
  template <typename Real> struct formatter<Utils::Poly<Real>> : ostream_formatter
  {
  };
  template <typename Real> struct formatter<Utils::Sturm<Real>> : ostream_formatter
  {
  };
}  // namespace fmt

#endif

//
// eof: Utils_Poly.hh
//
