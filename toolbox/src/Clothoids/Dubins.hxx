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
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

///
/// file: Dubins.hxx
///

namespace G2lib {

  /*\
   |   ____        _     _
   |  |  _ \ _   _| |__ (_)_ __  ___
   |  | | | | | | | '_ \| | '_ \/ __|
   |  | |_| | |_| | |_) | | | | \__ \
   |  |____/ \__,_|_.__/|_|_| |_|___/
  \*/

  //!
  //! Class to manage a circle arc
  //!
  class Dubins {

    CircleArc m_C1, m_C2, m_C3; //! Three arc solution of DUBINS problem

    //  {LSL, RSR, RSL, LSR, RLR, LRL}
    using DubinsType = enum class DubinsType : integer {
      LSL, RSR, RSL, LSR, RLR, LRL
    };

  public:

    //!
    //! Build an empty circle
    //!
    Dubins() = default;

    //!
    //! Build a copy of an existing Dubins problem.
    //!
    Dubins( Dubins const & s )
    { this->copy(s); }

    //!
    //! Construct a circle arc with the standard parameters.
    //!
    //! \param[in] x0     initial position x-coordinate
    //! \param[in] y0     initial position y-coordinate
    //! \param[in] theta0 initial angle
    //! \param[in] x1     final position x-coordinate
    //! \param[in] y1     final position y-coordinate
    //! \param[in] theta1 final angle
    //! \param[in] kmax   max curvature
    //!
    explicit
    Dubins(
      real_type x0,
      real_type y0,
      real_type theta0,
      real_type x1,
      real_type y1,
      real_type theta1,
      real_type k_max
    );

    //!
    //! Make a copy of an existing circle arc.
    //!
    void
    copy( Dubins const & d ) {
      m_C1 = d.m_C1;
      m_C2 = d.m_C2;
      m_C3 = d.m_C3;
    }
  };

}

///
/// eof: Dubins.hxx
///
