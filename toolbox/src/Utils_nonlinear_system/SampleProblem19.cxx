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

class SampleProblem19 : public NonlinearSystem
{
public:
  SampleProblem19() : NonlinearSystem( "Sample problem 19", "no doc", 2 ) {}

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = x( 0 ) * ( power2( x( 0 ) ) + power2( x( 1 ) ) );
    f( 1 ) = x( 1 ) * ( power2( x( 0 ) ) + power2( x( 1 ) ) );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    J.insert( 0, 0 ) = 3 * power2( x( 0 ) ) + power2( x( 1 ) );
    J.insert( 0, 1 ) = 2 * x( 0 ) * x( 1 );

    J.insert( 1, 0 ) = 2 * x( 0 ) * x( 1 );
    J.insert( 1, 1 ) = power2( x( 0 ) ) + 3 * power2( x( 1 ) );
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.setZero();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 3, 3;
  }
};
