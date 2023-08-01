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
  public:
    using DubinsType = enum class DubinsType : integer
    { LSL, RSR, LSR, RSL, LRL, RLR };

  private:
    DubinsType m_solution_type;

    CircleArc m_C0, m_C1, m_C2; //! Three arc solution of DUBINS problem

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
    //! Construct a Dubins solution
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
    ) {
      this->build( x0, y0, theta0, x1, y1, theta1, k_max );
    }

    //!
    //! Make a copy of an existing Dubins solution.
    //!
    void
    copy( Dubins const & d ) {
      m_C0 = d.m_C0;
      m_C1 = d.m_C1;
      m_C2 = d.m_C2;
      m_solution_type = d.m_solution_type;
    }

    //!
    //! Construct a Dubins solution
    //!
    //! \param[in] x0     initial position x-coordinate
    //! \param[in] y0     initial position y-coordinate
    //! \param[in] theta0 initial angle
    //! \param[in] x1     final position x-coordinate
    //! \param[in] y1     final position y-coordinate
    //! \param[in] theta1 final angle
    //! \param[in] kmax   max curvature
    //!
    bool
    build(
      real_type x0,
      real_type y0,
      real_type theta0,
      real_type x1,
      real_type y1,
      real_type theta1,
      real_type k_max
    );

    //!
    //! Return the first cicle of the Dubins solution
    //!
    CircleArc const & C0() const { return m_C0; }

    //!
    //! Return the second cicle of the Dubins solution
    //!
    CircleArc const & C1() const { return m_C1; }

    //!
    //! Return the third cicle of the Dubins solution
    //!
    CircleArc const & C2() const { return m_C2; }

    void
    get_solution( ClothoidList & CL ) const {
      CL.init();
      CL.reserve(3);
      CL.push_back( m_C0 );
      CL.push_back( m_C1 );
      CL.push_back( m_C2 );
    }

    real_type length() const { return m_C0.length()+m_C1.length()+m_C2.length(); }

    DubinsType solution_type() const { return m_solution_type; }

  };

  inline
  string
  to_string( Dubins::DubinsType n ) {
    string res = "";
    switch ( n ) {
    case Dubins::DubinsType::LSL: res = "LSL"; break;
    case Dubins::DubinsType::RSR: res = "RSR"; break;
    case Dubins::DubinsType::LSR: res = "LSR"; break;
    case Dubins::DubinsType::RSL: res = "RSL"; break;
    case Dubins::DubinsType::LRL: res = "LRL"; break;
    case Dubins::DubinsType::RLR: res = "RLR"; break;
    }
    return res;
  };
}

///
/// eof: Dubins.hxx
///
