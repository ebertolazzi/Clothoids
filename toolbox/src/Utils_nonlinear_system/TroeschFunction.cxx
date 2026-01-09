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

class TroeschFunction : public NonlinearSystem
{
  real_type const rho;
  real_type const h;

public:
  TroeschFunction( integer neq )
    : NonlinearSystem(
        "Troesch Function",
        "@inproceedings{Varadhan2009,\n"
        "  author={R. Varadhan and Paul D. Gilbert},\n"
        "  title={{BB:} An {R} Package for Solving a Large System of\n"
        "         Nonlinear Equations and for Optimizing a "
        "High-Dimensional\n"
        "         Nonlinear Objective Function},\n"
        "  year={2009}\n"
        "}\n",
        neq )
    , rho( 10 )
    , h( 1.0 / ( neq + 1 ) )
  {
    check_min_equations( n, 1 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type bf = rho * h * h;

    for ( integer i = 0; i < n; ++i ) f( i ) = 2 * x( i ) + bf * sinh( rho * x( i ) );

    f( 0 ) -= x( 1 );
    f( n - 1 ) -= x( n - 2 ) + 1;

    for ( integer i = 1; i < n - 1; ++i ) f( i ) -= x( i - 1 ) + x( i + 1 );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type bf = rho * rho * h * h;
    for ( integer i = 0; i < n; ++i ) J.insert( i, i ) = 2 + bf * cosh( rho * x( i ) );
    for ( integer i = 0; i < n - 1; ++i ) J.insert( i, i + 1 ) = -1;
    for ( integer i = 1; i < n; ++i ) J.insert( i, i - 1 ) = -1;
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    for ( integer i{ 0 }; i < n; ++i ) x0( i ) = ( ( 123 * i ) % 1001 ) / 1000.0;
  }
};
