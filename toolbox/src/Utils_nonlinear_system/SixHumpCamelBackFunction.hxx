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

class SixHumpCamelBackFunction : public NonlinearSystem
{
public:
  SixHumpCamelBackFunction()
    : NonlinearSystem(
        "Six Hump Camel Back function",
        "Molga M. and Smutnicki C. (2005).\n"
        "Test functions for optimization needs,\n"
        "http://www.zsd.ict.pwr.wroc.pl/files/docs/functions\n",
        2 )
  {
  }

  virtual void evaluate( Vector const & x_in, Vector & f ) const override
  {
    real_type x   = x_in[0];
    real_type y   = x_in[1];
    real_type t2  = x * x;
    real_type t5  = t2 * t2;
    real_type t10 = y * y;
    f( 0 )        = 8.0 * x - 0.84E1 * t2 * x + 2.0 * t5 * x + y;
    f( 1 )        = x - 8.0 * y + 16.0 * t10 * y;
  }

  virtual void jacobian( Vector const & x_in, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type x      = x_in[0];
    real_type y      = x_in[1];
    real_type t1     = x * x;
    real_type t3     = t1 * t1;
    real_type t6     = y * y;
    J.insert( 0, 0 ) = 8.0 - 0.252E2 * t1 + 10.0 * t3;
    J.insert( 0, 1 ) = J.insert( 1, 0 ) = 1.0;
    J.insert( 1, 1 )                    = -8.0 + 48.0 * t6;
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << -0.08984, 0.71266;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << -5, 5;
  }
};
