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

class ScalarProblem : public NonlinearSystem
{
public:
  ScalarProblem() : NonlinearSystem( "Scalar problem f(x) = x * ( x - 5 )**2", "no doc", 1 ) {}

  virtual void evaluate( Vector const & x, Vector & f ) const override { f( 0 ) = x( 0 ) * power2( x( 0 ) - 5 ); }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    J.insert( 0, 0 ) = ( 3 * x( 0 ) - 5 ) * ( x( 0 ) - 5 );
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 2 );
    auto & x0{ x_vec[0] };
    auto & x1{ x_vec[1] };
    x0.resize( n );
    x1.resize( n );
    x0( 0 ) = 0;
    x1( 0 ) = 5;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0( 0 ) = 1;
  }
};
