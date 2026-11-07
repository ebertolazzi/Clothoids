/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Pipal project is distributed under the MIT License.                                       *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef INCLUDE_PIPAL_INPUT_HXX
#define INCLUDE_PIPAL_INPUT_HXX

namespace Pipal
{
  /**
   * \brief Build the input data for the solver.
   *
   * \param[in] name Problem name.
   * \param[in] x0 Initial guess for the solution.
   * \param[in] bl Lower bounds for the variables.
   * \param[in] bu Upper bounds for the variables.
   * \param[in] cl Lower bounds for the constraints.
   * \param[in] cu Upper bounds for the constraints.
   */
  template <typename Real> void Solver<Real>::buildInput(
    std::string const &  name,
    Vector<Real> const & x0,
    Vector<Real> const & bl,
    Vector<Real> const & bu,
    Vector<Real> const & cl,
    Vector<Real> const & cu )
  {
    // Create alias for easier access
    Parameter<Real> & p{ this->m_parameter };
    Input<Real> &     i{ this->m_input };

    // Set problem identity
    i.name = name;

    // Set number of original formulation variables
    i.n0 = detail::checked_integer( x0.size(), "Pipal problem dimension" );

    // Find indices sets
    static constexpr Real EPSILON{ std::numeric_limits<Real>::epsilon() };
    Mask const            cond_bl( bl.array() <= -p.rhs_bnd );
    Mask const            cond_bu( bu.array() >= p.rhs_bnd );
    Mask const            cond_cl( cl.array() <= -p.rhs_bnd );
    Mask const            cond_cu( cu.array() >= p.rhs_bnd );
    Mask const            cond_beq( ( bl.array() - bu.array() ).abs() <= EPSILON );
    Mask const            cond_ceq( ( cl.array() - cu.array() ).abs() <= EPSILON );

    // Now the index sets
    i.I1 = find( cond_bl && cond_bu );
    i.I2 = find( cond_beq );
    i.I3 = find( !cond_bl && cond_bu );
    i.I4 = find( cond_bl && !cond_bu );
    i.I5 = find( !cond_bl && !cond_bu && !cond_beq );
    i.I6 = find( cond_ceq );
    i.I7 = find( !cond_cl && cond_cu );
    i.I8 = find( cond_cl && !cond_cu );
    i.I9 = find( !cond_cl && !cond_cu && !cond_ceq );

    // Set right-hand side values
    i.b2 = bl( i.I2 );
    i.l3 = bl( i.I3 );
    i.u4 = bu( i.I4 );
    i.l5 = bl( i.I5 );
    i.u5 = bu( i.I5 );
    i.b6 = cl( i.I6 );
    i.l7 = cl( i.I7 );
    i.u8 = cu( i.I8 );
    i.l9 = cl( i.I9 );
    i.u9 = cu( i.I9 );

    // Set sizes of indices sets
    i.n1 = detail::checked_integer( i.I1.size(), "Pipal variable count" );
    i.n2 = detail::checked_integer( i.I2.size(), "Pipal variable count" );
    i.n3 = detail::checked_integer( i.I3.size(), "Pipal variable count" );
    i.n4 = detail::checked_integer( i.I4.size(), "Pipal variable count" );
    i.n5 = detail::checked_integer( i.I5.size(), "Pipal variable count" );
    i.n6 = detail::checked_integer( i.I6.size(), "Pipal constraint count" );
    i.n7 = detail::checked_integer( i.I7.size(), "Pipal constraint count" );
    i.n8 = detail::checked_integer( i.I8.size(), "Pipal constraint count" );
    i.n9 = detail::checked_integer( i.I9.size(), "Pipal constraint count" );

    // Initialize number of invalid bounds
    i.vi = 0;

    // Count invalid bounds
    if ( i.n2 > 0 )
    {
      detail::add_count( i.vi, ( i.b2.array() <= -p.rhs_bnd ).count(), "Pipal invalid-bound count" );
      detail::add_count( i.vi, ( i.b2.array() >= p.rhs_bnd ).count(), "Pipal invalid-bound count" );
    }
    if ( i.n3 > 0 ) { detail::add_count( i.vi, ( i.l3.array() >= p.rhs_bnd ).count(), "Pipal invalid-bound count" ); }
    if ( i.n4 > 0 ) { detail::add_count( i.vi, ( i.u4.array() <= -p.rhs_bnd ).count(), "Pipal invalid-bound count" ); }
    if ( i.n5 > 0 )
    {
      detail::add_count( i.vi, ( i.l5.array() >= p.rhs_bnd ).count(), "Pipal invalid-bound count" );
      detail::add_count( i.vi, ( i.u5.array() <= -p.rhs_bnd ).count(), "Pipal invalid-bound count" );
      detail::add_count( i.vi, ( i.l5.array() > i.u5.array() ).count(), "Pipal invalid-bound count" );
    }
    if ( i.n6 > 0 )
    {
      detail::add_count( i.vi, ( i.b6.array() <= -p.rhs_bnd ).count(), "Pipal invalid-bound count" );
      detail::add_count( i.vi, ( i.b6.array() >= p.rhs_bnd ).count(), "Pipal invalid-bound count" );
    }
    if ( i.n7 > 0 ) { detail::add_count( i.vi, ( i.l7.array() >= p.rhs_bnd ).count(), "Pipal invalid-bound count" ); }
    if ( i.n8 > 0 ) { detail::add_count( i.vi, ( i.u8.array() <= -p.rhs_bnd ).count(), "Pipal invalid-bound count" ); }
    if ( i.n9 > 0 )
    {
      detail::add_count( i.vi, ( i.l9.array() >= p.rhs_bnd ).count(), "Pipal invalid-bound count" );
      detail::add_count( i.vi, ( i.u9.array() <= -p.rhs_bnd ).count(), "Pipal invalid-bound count" );
      detail::add_count( i.vi, ( i.l9.array() > i.u9.array() ).count(), "Pipal invalid-bound count" );
    }

    // Set number of variables and constraints
    i.nV = detail::checked_integer(
      static_cast<Eigen::Index>( i.n1 ) + i.n3 + i.n4 + i.n5, "Pipal internal variable count" );
    i.nI = detail::checked_integer(
      static_cast<Eigen::Index>( i.n3 ) + i.n4 + 2 * static_cast<Eigen::Index>( i.n5 ) + i.n7 + i.n8 +
        2 * static_cast<Eigen::Index>( i.n9 ),
      "Pipal internal constraint count" );
    i.nE = i.n6;

    // Set size of primal-dual matrix
    i.nA = detail::checked_integer(
      static_cast<Eigen::Index>( i.nV ) + 3 * static_cast<Eigen::Index>( i.nE ) + 3 * static_cast<Eigen::Index>( i.nI ),
      "Pipal primal-dual matrix size" );

    // Set initial point
    i.x0.resize( i.nV );
    i.x0 << x0( i.I1 ), x0( i.I3 ), x0( i.I4 ), x0( i.I5 );
  }

}  // namespace Pipal

#endif  // INCLUDE_PIPAL_INPUT_HXX
