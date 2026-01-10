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

class TridimensionalValley : public NonlinearSystem
{
  real_type const c1;
  real_type const c2;

public:
  TridimensionalValley()
    : NonlinearSystem( "Tridimensional valley.", "no doc", 3 ), c1( 1.003344481605351 ), c2( -3.344481605351171E-3 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type bf = ( c2 * power2( x( 0 ) ) + c1 ) * x( 0 );
    f( 0 )       = bf * exp( -power2( x( 0 ) ) / 100 ) - 1;
    f( 1 )       = 10 * ( sin( x( 0 ) ) - x( 1 ) );
    f( 2 )       = 10 * ( cos( x( 0 ) ) - x( 2 ) );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type t1     = x( 0 ) * x( 0 );
    J.insert( 0, 0 ) = ( c1 - ( t1 / 50.0 ) * ( c1 + c2 * ( t1 - 150.0 ) ) ) * exp( -t1 / 100.0 );
    J.insert( 1, 0 ) = 10 * cos( x( 0 ) );
    J.insert( 1, 1 ) = -10;
    J.insert( 2, 0 ) = -10 * sin( x( 0 ) );
    J.insert( 2, 2 ) = -10;
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 1.0103301175891008618821430258424435903873866121054053,
      0.847007375051043571769939585744456415641478463861070185,
      0.531581138312055623979884869864864195697816223034820704;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << -4, 1, 2;
  }
};
