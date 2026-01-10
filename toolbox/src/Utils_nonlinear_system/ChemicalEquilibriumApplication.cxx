/*\
 |
 |  Author:
 |    Enrico Bertolazzi
 |    University of Trento
 |    Department of Industrial Engineering
 |    Via Sommarive 9, I-38123, Povo, Trento, Italy
 |    email: enrico.bertolazzi@unitn.it
\*/

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class ChemicalEquilibriumApplication : public NonlinearSystem
{
  real_type const R;
  real_type const R5;
  real_type const R6;
  real_type const R7;
  real_type const R8;
  real_type const R9;
  real_type const R10;

public:
  ChemicalEquilibriumApplication()
    : NonlinearSystem(
        "Chemical Equilibrium Application",
        "@article{Meintjes:1990,\n"
        "  author  = {Meintjes, Keith and Morgan, Alexander P.},\n"
        "  title   = {Chemical Equilibrium Systems As Numerical Test "
        "Problems},\n"
        "  journal = {ACM Trans. Math. Softw.},\n"
        "  year    = {1990},\n"
        "  volume  = {16},\n"
        "  number  = {2},\n"
        "  pages   = {143--151},\n"
        "  doi     = {10.1145/78928.78930},\n"
        "}\n\n"
        "@article{Hentenryck:1997,\n"
        "  author  = {Van Hentenryck, P. and McAllester, D. and Kapur, "
        "D.},\n"
        "  title   = {Solving Polynomial Systems Using a Branch and Prune "
        "Approach},\n"
        "  journal = {SIAM Journal on Numerical Analysis},\n"
        "  year    = {1997},\n"
        "  volume  = {34},\n"
        "  number  = {2},\n"
        "  pages   = {797-827},\n"
        "  doi = {10.1137/S0036142995281504}\n"
        "}\n",
        5 )
    , R( 10 )
    , R5( 0.193 )
    , R6( 0.002597 / sqrt( 40.0 ) )
    , R7( 0.003448 / sqrt( 40.0 ) )
    , R8( 0.00001799 / sqrt( 40.0 ) )
    , R9( 0.0002155 / sqrt( 40.0 ) )
    , R10( 0.00003846 / sqrt( 40.0 ) )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = x( 0 ) * x( 1 ) + x( 0 ) - 3 * x( 4 );
    f( 1 ) = 2 * x( 0 ) * x( 1 ) + x( 0 ) + x( 1 ) * x( 2 ) * x( 2 ) + R8 * x( 1 ) - R * x( 4 ) +
             2 * R10 * x( 1 ) * x( 1 ) + R7 * x( 1 ) * x( 2 ) + R9 * x( 1 ) * x( 3 );
    f( 2 ) = 2 * ( x( 1 ) + R5 ) * x( 2 ) * x( 2 ) - 8 * x( 4 ) + R6 * x( 2 ) + R7 * x( 1 ) * x( 2 );
    f( 3 ) = ( R9 * x( 1 ) + 2 * x( 3 ) ) * x( 3 ) - 4 * R * x( 4 );
    f( 4 ) = x( 0 ) * ( x( 1 ) + 1 ) + ( R8 + R10 * x( 1 ) + R7 * x( 2 ) + R9 * x( 3 ) ) * x( 1 ) +
             ( R5 + x( 1 ) ) * x( 2 ) * x( 2 ) + x( 3 ) * x( 3 ) - 1 + R6 * x( 2 );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    J.insert( 0, 0 ) = 1 + x( 1 );
    J.insert( 0, 1 ) = x( 0 );
    J.insert( 0, 4 ) = -3;

    J.insert( 1, 0 ) = 2 * x( 1 ) + 1;
    J.insert( 1, 1 ) = 4 * R10 * x( 1 ) + R7 * x( 2 ) + R9 * x( 3 ) + x( 2 ) * x( 2 ) + R8 + 2 * x( 0 );
    J.insert( 1, 2 ) = 2 * x( 1 ) * x( 2 ) + R7 * x( 1 );
    J.insert( 1, 3 ) = R9 * x( 1 );
    J.insert( 1, 4 ) = -R;

    J.insert( 2, 1 ) = 2 * x( 2 ) * x( 2 ) + R7 * x( 2 );
    J.insert( 2, 2 ) = 4 * ( x( 1 ) + R5 ) * x( 2 ) + R6 + R7 * x( 1 );
    J.insert( 2, 4 ) = -8;

    J.insert( 3, 1 ) = R9 * x( 3 );
    J.insert( 3, 3 ) = R9 * x( 1 ) + 4 * x( 3 );
    J.insert( 3, 4 ) = -4 * R;

    J.insert( 4, 0 ) = x( 1 ) + 1;
    J.insert( 4, 1 ) = x( 0 ) + R8 + 2 * R10 * x( 1 ) + R7 * x( 2 ) + R9 * x( 3 ) + x( 2 ) * x( 2 );
    J.insert( 4, 2 ) = R7 * x( 1 ) + 2 * ( R5 + x( 1 ) ) * x( 2 ) + R6;
    J.insert( 4, 3 ) = R9 * x( 1 ) + 2 * x( 3 );

    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 1 );
  }
};
