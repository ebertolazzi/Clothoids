
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

class BertolazziRootPlusSquare : public NonlinearSystem
{
  real_type const epsilon;
  real_type const delta;
  real_type       xmin;

public:
  BertolazziRootPlusSquare()
    : NonlinearSystem( "Bertolazzi: root+square.", "no doc", 2 ), epsilon( 1E-6 ), delta( 1E-8 )
  {
    xmin = epsilon * ( 1.0 / exp( 1.0 ) - 1.0 );
  }

  virtual void evaluate( Vector const & x_in, Vector & f ) const override
  {
    real_type x = x_in( 0 );
    real_type y = x_in( 1 );
    f( 0 )      = x > xmin ? log1p( log1p( x / epsilon ) ) : nan( "nan" );
    f( 1 )      = delta * y + power2( y );
  }

  virtual void jacobian( Vector const & x_in, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type x      = x_in( 0 );
    real_type y      = x_in( 1 );
    J.insert( 0, 0 ) = x > xmin ? 1 / ( ( x + epsilon ) * ( 1 + log1p( x / epsilon ) ) ) : nan( "nan" );
    J.insert( 1, 1 ) = delta + 2 * y;
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
    x0 << 1000, 1000;
  }

  virtual void check_if_admissible( Vector const & x ) const override
  {
    UTILS_ASSERT( x( 0 ) >= xmin, "BertolazziRootPlusSquare: x = {} must be >= {}", x( 0 ), xmin );
  }

  virtual void bounding_box( Vector & L, Vector & U ) const override
  {
    L( 0 ) = xmin;
    L( 1 ) = -real_max;
    U( 0 ) = U( 1 ) = real_max;
  }
};

/*
 * Reference:
 *  E.Bertolazzi
 */

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class BertolazziAtanPlusQuadratic : public NonlinearSystem
{
  real_type const epsilon;

public:
  BertolazziAtanPlusQuadratic() : NonlinearSystem( "Bertolazzi: atan+quadratic", "no doc", 2 ), epsilon( 1E-9 ) {}

  virtual void evaluate( Vector const & x_in, Vector & f ) const override
  {
    real_type x = x_in( 0 );
    real_type y = x_in( 1 );
    f( 0 )      = atan( x / epsilon );
    f( 1 )      = epsilon * y + x * y;
  }

  virtual void jacobian( Vector const & x_in, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type x      = x_in( 0 );
    real_type y      = x_in( 1 );
    real_type xe     = x / epsilon;
    J.insert( 0, 0 ) = 1 / ( epsilon + x * xe );
    J.insert( 1, 0 ) = y;
    J.insert( 1, 1 ) = epsilon + x;
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
    x0 << 1, 1;
  }

  virtual void check_if_admissible( Vector const & x ) const override
  {
    for ( integer i = 0; i < n; ++i )
      UTILS_ASSERT( std::abs( x( i ) ) < 10, "x[{}] = {} out of range [-10,10]", i, x( i ) );
  }

  virtual void bounding_box( Vector & L, Vector & U ) const override
  {
    L.fill( -10 );
    U.fill( 10 );
  }
};

/*
 * Reference:
 *  E.Bertolazzi
 *
 * f = x*(1+x^2)*exp(-x)+exp(x/10)-1
 * g = x*exp(-x)+exp(x/10)-(3*exp(-3)+exp(3/10));
 *
 * f-g
 * f+g
 */

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class BertolazziHard : public NonlinearSystem
{
public:
  BertolazziHard() : NonlinearSystem( "Bertolazzi Hard.", "no doc", 2 ) {}

  virtual void evaluate( Vector const & x_in, Vector & f ) const override
  {
    real_type x  = x_in( 0 );
    real_type y  = x_in( 1 );
    real_type t1 = exp( x / 10 );
    real_type t2 = exp( y / 10 );
    real_type t3 = exp( 3.0 / 10.0 );
    real_type t4 = x * ( x * x + 1 ) * exp( -x );
    real_type t5 = y * exp( -y );
    real_type t6 = 3 * exp( -3 );
    f( 0 )       = t6 - 1 + t4 - t5 + t1 - t2 + t3;
    f( 1 )       = -t6 - 1 + t4 + t5 + t1 + t2 - t3;
  }

  virtual void jacobian( Vector const & x_in, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type x      = x_in( 0 );
    real_type y      = x_in( 1 );
    real_type t1     = x * x;
    t1               = exp( x / 10 ) / 10 + exp( -x ) * ( 3 * t1 + 1 - x * ( 1 + t1 ) );
    real_type t2     = -exp( y / 10 ) / 10 + exp( -y ) * ( y - 1 );
    J.insert( 0, 0 ) = t1;
    J.insert( 0, 1 ) = t2;
    J.insert( 1, 0 ) = t1;
    J.insert( 1, 1 ) = -t2;
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 0, 3;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 10, -1;
  }
};

/*
 * Reference:
 *  E.Bertolazzi
 *
 * f = x*(1+x^2)*exp(-x)+exp(x/10)-1
 */

class BertolazziSingleEQ : public NonlinearSystem
{
public:
  BertolazziSingleEQ() : NonlinearSystem( "Bertolazzi Single EQ.", "no doc", 1 ) {}

  virtual void evaluate( Vector const & x_in, Vector & f ) const override
  {
    real_type x = x_in( 0 );
    f( 0 )      = exp( -x ) * ( x * x + 1 ) * x + exp( x / 10 ) - 1;
  }

  virtual void jacobian( Vector const & x_in, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type x      = x_in( 0 );
    real_type t1     = x * x;
    real_type t2     = exp( -x );
    real_type t9     = exp( x / 10 );
    J.insert( 0, 0 ) = 3 * t2 * t1 + t2 - t1 * x * t2 - t2 * x + t9 / 10;
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 0;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 10;
  }
};
