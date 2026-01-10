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

/*
Numerical Performance of Abs Codes for Nonlinear Systems of Equations
E. Bodon (1), A. Del Popolo (1, 2, 3), L. Luksan (4), E. Spedicato (1)
2001
https://arxiv.org/abs/math/0106029
*/

class SchubertBroydenFunction : public NonlinearSystem
{
public:
  SchubertBroydenFunction( integer neq ) : NonlinearSystem( "Schubert Broyden function", "no doc", neq ) {}

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 )     = ( 3 - x( 0 ) ) * x( 0 ) + 1 - 2 * x( 1 );
    f( n - 1 ) = ( 3 - x( n - 1 ) ) * x( n - 1 ) + 1 - 2 * x( n - 2 );
    for ( integer i = 1; i < n - 1; ++i ) f( i ) = ( 3 - x( i ) ) * x( i ) - x( i - 1 ) - 2 * x( i + 1 );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    // Prima equazione (i=0): (3-x(0))*x(0) + 1 - 2*x(1)
    J.insert( 0, 0 ) = 3 - 2 * x( 0 );  // derivata di (3-x(0))*x(0) = 3x(0) - x(0)^2
    J.insert( 0, 1 ) = -2;              // derivata di -2*x(1)

    // Ultima equazione (i=n-1): (3-x(n-1))*x(n-1) + 1 - 2*x(n-2)
    J.insert( n - 1, n - 1 ) = 3 - 2 * x( n - 1 );  // derivata di (3-x(n-1))*x(n-1)
    J.insert( n - 1, n - 2 ) = -2;                  // derivata di -2*x(n-2)

    // Equazioni intermedie (i da 1 a n-2): (3-x(i))*x(i) - x(i-1) - 2*x(i+1)
    for ( integer i = 1; i < n - 1; ++i )
    {
      J.insert( i, i - 1 ) = -1;              // derivata di -x(i-1)
      J.insert( i, i )     = 3 - 2 * x( i );  // derivata di (3-x(i))*x(i)
      J.insert( i, i + 1 ) = -2;              // derivata di -2*x(i+1)
    }

    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( -1 );
  }
};
