/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2025                                                      |
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

#ifndef UTILS_NONLINEAR_SYSTEM_HH
#define UTILS_NONLINEAR_SYSTEM_HH

#include "Utils.hh"
#include "Utils_eigen.hh"
#include "Utils_fmt.hh"

//! namespace for nonlinear systems and nonlinearsolver
namespace Utils
{

  using namespace ::std;

  void printNotConverged( std::basic_ostream<char> & stream );
  void printConverged( std::basic_ostream<char> & stream );

  /*
  //   _   _             _ _                       ____            _
  //  | \ | | ___  _ __ | (_)_ __   ___  __ _ _ __/ ___| _   _ ___| |_ ___ _ __
  ___
  //  |  \| |/ _ \| '_ \| | | '_ \ / _ \/ _` | '__\___ \| | | / __| __/ _ \ '_ `
  _ \
  //  | |\  | (_) | | | | | | | | |  __/ (_| | |   ___) | |_| \__ \ ||  __/ | |
  | | |
  //  |_| \_|\___/|_| |_|_|_|_| |_|\___|\__,_|_|  |____/ \__, |___/\__\___|_|
  |_| |_|
  //                                                     |___/
  */

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! base class for nonlinear system
  class NonlinearSystem
  {
  public:
    using integer      = Eigen::Index;
    using real_type    = double;
    using Vector       = Eigen::Matrix<real_type, Eigen::Dynamic, 1>;
    using SparseMatrix = Eigen::SparseMatrix<real_type>;

    NonlinearSystem( NonlinearSystem const & )                   = delete;
    NonlinearSystem const & operator=( NonlinearSystem const & ) = delete;

  private:
    string const m_title;
    string const m_bibtex;

  protected:
    void
    check_min_equations( integer i, integer i_min ) const
    {
      UTILS_ASSERT( i >= i_min, "check_min_equations:: i = {} < {}", i, i_min );
    }

    void
    check_even( integer i, integer i_min ) const
    {
      UTILS_ASSERT( ( i % 2 ) == 0 && i >= i_min, "check_even:: odd index i = {}" );
    }

    void
    check_odd( integer i, integer i_min ) const
    {
      UTILS_ASSERT( ( i % 2 ) != 0 && i >= i_min, "check_odd:: odd index i = {}" );
    }

    void
    check_three( integer i, integer i_min ) const
    {
      UTILS_ASSERT( ( i % 3 ) == 0 && i >= i_min, "check_three:: index i = {}" );
    }

    void
    check_four( integer i, integer i_min ) const
    {
      UTILS_ASSERT( ( i % 4 ) == 0 && i >= i_min, "check_four:: index i = {}" );
    }

    real_type real_max{ std::numeric_limits<real_type>::max() };

    integer n;

  public:
    NonlinearSystem( string const & title, string const & bibtex, integer dim )
      : m_title( fmt::format( "{} neq={}", title, dim ) ), m_bibtex( bibtex ), n( dim )
    {
    }

    virtual ~NonlinearSystem() {}

    //!
    //! \name Virtual functions common for all nonlinear systems
    //!
    //! @{
    //!
    //! The number of equations of the nonlinear system
    //!
    integer
    num_equations( void ) const
    {
      return n;
    }

    //!
    //! Check if `x` is an acceptable point
    //!
    virtual void check_if_admissible( Vector const & ) const {};

    //!
    //! Compute \f$ F(x) \f$
    //!
    virtual void evaluate( Vector const & x, Vector & f ) const = 0;

    //!
    //! Evaluate \f$ J(x) = \dfrac{}{} F(x) \f$
    //!
    virtual void jacobian( Vector const & x, SparseMatrix & jac ) const = 0;

    //!
    //! Box where to search the solution
    //!
    virtual void
    bounding_box( Vector & L, Vector & U ) const
    {
      L.fill( -real_max );
      U.fill( real_max );
    }

    //!
    //! Get initial guesses for the search of solution
    //!
    virtual void initial_points( vector<Vector> & x_vec ) const = 0;

    virtual void
    exact_solution( vector<Vector> & x_vec ) const
    {
      x_vec.clear();
    }

    //!
    //! @}
    //!

    //! Bibliography of the test
    string const &
    bibtex() const
    {
      return m_bibtex;
    }

    //! The name of the test
    string const &
    title() const
    {
      return m_title;
    }
  };

  extern std::vector<NonlinearSystem *> nonlinear_system_tests;
  extern std::map<string, unsigned>     nonlinear_system_tests_map;
  extern void                           init_nonlinear_system_tests();

}  // namespace Utils

#endif
