/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                      |
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

#include "Clothoids.hh"
#include "Clothoids_fmt.hh"

namespace G2lib
{

#ifdef AUTODIFF_SUPPORT

  autodiff::dual1st ClothoidData::X( autodiff::dual1st const & s ) const
  {
    real_type C, S;
    real_type ss = s.val;
    GeneralizedFresnelCS( m_dk * ( ss * ss ), m_kappa0 * ss, m_theta0, C, S );
    autodiff::dual1st res;
    res.val  = m_x0 + ss * C;
    res.grad = s.grad * cos( theta( ss ) );
    return res;
  }

  autodiff::dual1st ClothoidData::Y( autodiff::dual1st const & s ) const
  {
    real_type C, S;
    real_type ss = s.val;
    GeneralizedFresnelCS( m_dk * ( ss * ss ), m_kappa0 * ss, m_theta0, C, S );
    autodiff::dual1st res;
    res.val  = m_y0 + ss * S;
    res.grad = s.grad * sin( theta( ss ) );  // Added missing semicolon
    return res;
  }

  autodiff::dual2nd ClothoidData::X( autodiff::dual2nd const & s ) const
  {
    real_type C, S;
    real_type ss = s.val.val;
    GeneralizedFresnelCS( m_dk * ( ss * ss ), m_kappa0 * ss, m_theta0, C, S );

    autodiff::dual2nd res;

    // Value
    real_type V = m_x0 + ss * C;

    // First and second derivatives
    real_type const sdk       = ss * m_dk;
    real_type const theta_val = m_theta0 + ss * ( m_kappa0 + 0.5 * sdk );
    real_type const theta_D   = m_kappa0 + sdk;
    real_type       D         = cos( theta_val );
    real_type       DD        = -sin( theta_val ) * theta_D;

    // Set value
    res.val.val = V;

    // Set first derivative (gradient of value)
    // f'(x) = cos(theta(s)) * s'
    res.val.grad = s.val.grad * D;

    // Set gradient (first derivative) and its derivative
    // grad.val = f'(x) = cos(theta(s)) * s'
    res.grad.val = D * s.grad.val;

    // grad.grad = f''(x) = [cos(theta(s)) * s']'
    //           = -sin(theta(s)) * theta'(s) * (s')^2 + cos(theta(s)) * s''
    // For independent variable s, s' = 1, s'' = 0
    res.grad.grad = DD * s.val.grad * s.grad.val + D * s.grad.grad;

    return res;
  }

  autodiff::dual2nd ClothoidData::Y( autodiff::dual2nd const & s ) const
  {
    real_type C, S;
    real_type ss = s.val.val;
    GeneralizedFresnelCS( m_dk * ( ss * ss ), m_kappa0 * ss, m_theta0, C, S );

    autodiff::dual2nd res;

    // Value
    real_type V = m_y0 + ss * S;

    // First and second derivatives
    real_type const sdk       = ss * m_dk;  // Fixed typo: was `s * m_dk`
    real_type const theta_val = m_theta0 + ss * ( m_kappa0 + 0.5 * sdk );
    real_type const theta_D   = m_kappa0 + sdk;
    real_type       D         = sin( theta_val );
    real_type       DD        = cos( theta_val ) * theta_D;

    // Set value
    res.val.val = V;

    // Set first derivative (gradient of value)
    // f'(x) = sin(theta(s)) * s'
    res.val.grad = s.val.grad * D;

    // Set gradient (first derivative) and its derivative
    // grad.val = f'(x) = sin(theta(s)) * s'
    res.grad.val = D * s.grad.val;

    // grad.grad = f''(x) = [sin(theta(s)) * s']'
    //           = cos(theta(s)) * theta'(s) * (s')^2 + sin(theta(s)) * s''
    // For independent variable s, s' = 1, s'' = 0
    res.grad.grad = DD * s.val.grad * s.grad.val + D * s.grad.grad;

    return res;
  }

  // Tangent and normal vectors for clothoid curve
  // theta(s) = theta0 + kappa0 * s + 0.5 * dk * s^2
  // Derivative: theta'(s) = kappa0 + dk * s

  // Tangent x-component (cos(theta))
  autodiff::dual1st ClothoidData::tg_x( autodiff::dual1st const & s ) const
  {
    real_type ss        = s.val;
    real_type theta_val = m_theta0 + m_kappa0 * ss + 0.5 * m_dk * ss * ss;
    real_type theta_der = m_kappa0 + m_dk * ss;

    autodiff::dual1st res;
    res.val  = cos( theta_val );
    res.grad = -sin( theta_val ) * theta_der * s.grad;
    return res;
  }

  autodiff::dual2nd ClothoidData::tg_x( autodiff::dual2nd const & s ) const
  {
    real_type ss = s.val.val;
    real_type s1 = s.val.grad;
    real_type s2 = s.grad.grad;

    real_type theta_val  = m_theta0 + m_kappa0 * ss + 0.5 * m_dk * ss * ss;
    real_type theta_der1 = m_kappa0 + m_dk * ss;

    // First derivative: f'(s) = -sin(theta) * theta'
    // Second derivative: f''(s) = -cos(theta) * (theta')^2 - sin(theta) * theta''
    // where theta'' = dk
    real_type f_val  = cos( theta_val );
    real_type f_der1 = -sin( theta_val ) * theta_der1;
    real_type f_der2 = -cos( theta_val ) * theta_der1 * theta_der1 - sin( theta_val ) * m_dk;

    autodiff::dual2nd res;

    // Value
    res.val.val  = f_val;
    res.val.grad = f_der1 * s1;

    // Gradient (first derivative)
    res.grad.val = f_der1 * s.grad.val;

    // Second derivative: f''(s) * (s1)^2 + f'(s) * s2
    // For independent variable s: s1 = 1, s2 = 0
    res.grad.grad = f_der2 * s1 * s.grad.val + f_der1 * s2;

    return res;
  }

  // Tangent y-component (sin(theta))
  autodiff::dual1st ClothoidData::tg_y( autodiff::dual1st const & s ) const
  {
    real_type ss        = s.val;
    real_type theta_val = m_theta0 + m_kappa0 * ss + 0.5 * m_dk * ss * ss;
    real_type theta_der = m_kappa0 + m_dk * ss;

    autodiff::dual1st res;
    res.val  = sin( theta_val );
    res.grad = cos( theta_val ) * theta_der * s.grad;
    return res;
  }

  autodiff::dual2nd ClothoidData::tg_y( autodiff::dual2nd const & s ) const
  {
    real_type ss = s.val.val;
    real_type s1 = s.val.grad;
    real_type s2 = s.grad.grad;

    real_type theta_val  = m_theta0 + m_kappa0 * ss + 0.5 * m_dk * ss * ss;
    real_type theta_der1 = m_kappa0 + m_dk * ss;

    // First derivative: f'(s) = cos(theta) * theta'
    // Second derivative: f''(s) = -sin(theta) * (theta')^2 + cos(theta) * theta''
    real_type f_val  = sin( theta_val );
    real_type f_der1 = cos( theta_val ) * theta_der1;
    real_type f_der2 = -sin( theta_val ) * theta_der1 * theta_der1 + cos( theta_val ) * m_dk;

    autodiff::dual2nd res;

    // Value
    res.val.val  = f_val;
    res.val.grad = f_der1 * s1;

    // Gradient
    res.grad.val = f_der1 * s.grad.val;

    // Second derivative
    res.grad.grad = f_der2 * s1 * s.grad.val + f_der1 * s2;

    return res;
  }

  // Normal vector components (ISO standard: right-handed system)
  // nor_x_ISO = -sin(theta), nor_y_ISO = cos(theta)
  autodiff::dual1st ClothoidData::nor_x_ISO( autodiff::dual1st const & s ) const
  {
    real_type ss        = s.val;
    real_type theta_val = m_theta0 + m_kappa0 * ss + 0.5 * m_dk * ss * ss;
    real_type theta_der = m_kappa0 + m_dk * ss;

    autodiff::dual1st res;
    res.val  = -sin( theta_val );
    res.grad = -cos( theta_val ) * theta_der * s.grad;
    return res;
  }

  autodiff::dual2nd ClothoidData::nor_x_ISO( autodiff::dual2nd const & s ) const
  {
    real_type ss = s.val.val;
    real_type s1 = s.val.grad;
    real_type s2 = s.grad.grad;

    real_type theta_val  = m_theta0 + m_kappa0 * ss + 0.5 * m_dk * ss * ss;
    real_type theta_der1 = m_kappa0 + m_dk * ss;

    // First derivative: f'(s) = -cos(theta) * theta'
    // Second derivative: f''(s) = sin(theta) * (theta')^2 - cos(theta) * theta''
    real_type f_val  = -sin( theta_val );
    real_type f_der1 = -cos( theta_val ) * theta_der1;
    real_type f_der2 = sin( theta_val ) * theta_der1 * theta_der1 - cos( theta_val ) * m_dk;

    autodiff::dual2nd res;

    res.val.val   = f_val;
    res.val.grad  = f_der1 * s1;
    res.grad.val  = f_der1 * s.grad.val;
    res.grad.grad = f_der2 * s1 * s.grad.val + f_der1 * s2;

    return res;
  }

  autodiff::dual1st ClothoidData::nor_y_ISO( autodiff::dual1st const & s ) const
  {
    real_type ss        = s.val;
    real_type theta_val = m_theta0 + m_kappa0 * ss + 0.5 * m_dk * ss * ss;
    real_type theta_der = m_kappa0 + m_dk * ss;

    autodiff::dual1st res;
    res.val  = cos( theta_val );
    res.grad = -sin( theta_val ) * theta_der * s.grad;
    return res;
  }

  autodiff::dual2nd ClothoidData::nor_y_ISO( autodiff::dual2nd const & s ) const
  {
    real_type ss = s.val.val;
    real_type s1 = s.val.grad;
    real_type s2 = s.grad.grad;

    real_type theta_val  = m_theta0 + m_kappa0 * ss + 0.5 * m_dk * ss * ss;
    real_type theta_der1 = m_kappa0 + m_dk * ss;

    // Same as tg_x but with sign change in second derivative
    real_type f_val  = cos( theta_val );
    real_type f_der1 = -sin( theta_val ) * theta_der1;
    real_type f_der2 = -cos( theta_val ) * theta_der1 * theta_der1 - sin( theta_val ) * m_dk;

    autodiff::dual2nd res;

    res.val.val   = f_val;
    res.val.grad  = f_der1 * s1;
    res.grad.val  = f_der1 * s.grad.val;
    res.grad.grad = f_der2 * s1 * s.grad.val + f_der1 * s2;

    return res;
  }

#endif

}  // namespace G2lib

// EOF: Fresnel.cc
