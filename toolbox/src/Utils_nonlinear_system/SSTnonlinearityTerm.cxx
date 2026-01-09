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

static inline string ini_msg_SSTnonlinearityTerm( int item )
{
  return fmt::format( "SST nonlinearity term, N.{}", item );
}

class SSTnonlinearityTerm : public NonlinearSystem
{
  integer const   idx;
  real_type const SCALE;
  real_type       sst1;

public:
  SSTnonlinearityTerm( integer idx_in )
    : NonlinearSystem(
        ini_msg_SSTnonlinearityTerm( idx_in ),
        "@article{Sincovec:1975,\n"
        "  author  = {Sincovec, Richard F. and Madsen, Niel K.},\n"
        "  title   = {Software for Nonlinear Partial Differential "
        "Equations},\n"
        "  journal = {ACM Trans. Math. Softw.},\n"
        "  volume  = {1},\n"
        "  number  = {3},\n"
        "  year    = {1975},\n"
        "  pages   = {232--260},\n"
        "  doi     = {10.1145/355644.355649}\n"
        "}\n",
        4 )
    , idx( idx_in )
    , SCALE( 1e7 )
  {
    UTILS_ASSERT( idx == 0 || idx == 1, "SSTnonlinearityTerm, idx = {} must be 0 or 1", idx_in );
    switch ( idx )
    {
      case 0: sst1 = 360; break;
      case 1: sst1 = 3250; break;
    }
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = 4E5 - 272.443800016 * x( 0 ) + 0.1e-3 * x( 1 ) + 0.7e-2 * x( 3 ) - 3.67E-16 * x( 0 ) * x( 1 ) -
             4.13E-12 * x( 0 ) * x( 3 );
    f( 1 ) = 272.4438 * x( 0 ) - 0.100016e-3 * x( 1 ) + 3.67E-16 * x( 0 ) * x( 1 ) - 3.57E-15 * x( 1 ) * x( 2 );
    f( 2 ) = -1.6E-8 * x( 2 ) + 0.7e-2 * x( 3 ) + 4.1283E-12 * x( 0 ) * x( 3 ) - 3.57E-15 * x( 1 ) * x( 2 ) + 800.0 +
             sst1;
    f( 3 ) = -0.7000016e-2 * x( 3 ) + 3.57E-15 * x( 1 ) * x( 2 ) - 4.1283E-12 * x( 0 ) * x( 3 ) + 800.0;
    for ( integer i = 0; i < n; ++i ) f( i ) /= SCALE;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    J.insert( 0, 0 ) = -272.443800016 - 3.67E-16 * x( 1 ) - 4.13E-12 * x( 3 );
    J.insert( 0, 1 ) = 0.0001 - 3.67E-16 * x( 0 );
    J.insert( 0, 2 ) = 0.0;
    J.insert( 0, 3 ) = 0.007 - 4.13E-12 * x( 0 );

    J.insert( 1, 0 ) = 272.4438 + 3.67E-16 * x( 1 );
    J.insert( 1, 1 ) = -0.100016e-3 + 3.67E-16 * x( 0 ) - 3.57E-15 * x( 2 );
    J.insert( 1, 2 ) = -3.57E-15 * x( 1 );
    J.insert( 1, 3 ) = 0.0;

    J.insert( 2, 0 ) = 4.1283E-12 * x( 3 );
    J.insert( 2, 1 ) = -3.57E-15 * x( 2 );
    J.insert( 2, 2 ) = -1.6E-8 - 3.57E-15 * x( 1 );
    J.insert( 2, 3 ) = 0.007 + 4.1283E-12 * x( 0 );

    J.insert( 3, 0 ) = -4.1283E-12 * x( 3 );
    J.insert( 3, 1 ) = 3.57E-15 * x( 2 );
    J.insert( 3, 2 ) = 3.57E-15 * x( 1 );
    J.insert( 3, 3 ) = -0.7000016e-2 - 4.1283E-12 * x( 0 );

    J.makeCompressed();

    J /= SCALE;
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 2 );
    auto & x0{ x_vec[0] };
    auto & x1{ x_vec[1] };
    x0.resize( n );
    x1.resize( n );
    x0 << 1.264224341494800168746220557617102243320e6, 8.501430437421119197727722658565357159210e11,
      8.547006912362613482744674003963465220779e10, 3.702993087637386517255325996036534779221e10;
    x1 << 1.167055115589686729437819269080948198165e6, 3.069755937686399353641443136774453022013e11,
      2.621168343680655428727664834155029399049e11, 4.100816563193445712723351658449706009508e10;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 1E9, 1E9, 1E13, 1E7;
  }

  virtual void check_if_admissible( Vector const & x ) const override
  {
    for ( integer i = 0; i < n; ++i ) UTILS_ASSERT( x( i ) > 0, "Bad range" );
  }

  virtual void bounding_box( Vector & L, Vector & U ) const override
  {
    U.fill( real_max );
    L.setZero();
  }
};
