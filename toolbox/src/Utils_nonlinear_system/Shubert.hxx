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

class Shubert : public NonlinearSystem
{
public:
  Shubert()
    : NonlinearSystem(
        "Shubert Function",
        "@book{brent2013,\n"
        "  author    = {Brent, R.P.},\n"
        "  title     = {Algorithms for Minimization Without Derivatives},\n"
        "  isbn      = {9780486143682},\n"
        "  series    = {Dover Books on Mathematics},\n"
        "  year      = {2013},\n"
        "  publisher = {Dover Publications}\n"
        "}\n",
        2 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = f( 1 ) = 0;

    real_type factor1 = 0.0;
    real_type df1dx1  = 0.0;
    for ( integer i = 0; i < 5; ++i )
    {
      real_type y = i + 1;
      factor1 += y * cos( ( y + 1.0 ) * x( 0 ) + y );
      df1dx1 -= y * ( y + 1.0 ) * sin( ( y + 1.0 ) * x( 0 ) + y );
    }

    real_type factor2 = 0.0;
    real_type df2dx2  = 0.0;
    for ( integer i = 0; i < 5; ++i )
    {
      real_type y = i;
      factor2 += y * cos( ( y + 1.0 ) * x( 1 ) + y );
      df2dx2 -= y * ( y + 1.0 ) * sin( ( y + 1.0 ) * x( 1 ) + y );
    }
    f( 0 ) = df1dx1 * factor2;
    f( 1 ) = factor1 * df2dx2;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    real_type factor1   = 0.0;
    real_type df1dx1    = 0.0;
    real_type factor1_D = 0.0;
    real_type df1dx1_D  = 0.0;
    for ( integer i = 0; i < 5; ++i )
    {
      real_type y  = i + 1;
      real_type y1 = y + 1.0;
      factor1 += y * cos( y1 * x( 0 ) + y );
      df1dx1 -= y * y1 * sin( y1 * x( 0 ) + y );

      factor1_D -= y * sin( y1 * x( 0 ) + y ) * y1;
      df1dx1_D -= y * power2( y1 ) * cos( y1 * x( 0 ) + y );
    }

    real_type factor2   = 0.0;
    real_type df2dx2    = 0.0;
    real_type factor2_D = 0.0;
    real_type df2dx2_D  = 0.0;
    for ( integer i = 0; i < 5; ++i )
    {
      real_type y  = i;
      real_type y1 = y + 1.0;
      factor2 += y * cos( y1 * x( 1 ) + y );
      df2dx2 -= y * y1 * sin( y1 * x( 1 ) + y );
      factor2_D -= y * sin( y1 * x( 1 ) + y ) * y1;
      df2dx2_D -= y * power2( y1 ) * cos( y1 * x( 1 ) + y );
    }
    J.resize( n, n );
    J.setZero();
    J.insert( 0, 0 ) = df1dx1_D * factor2;
    J.insert( 0, 1 ) = df1dx1 * factor2_D;
    J.insert( 1, 0 ) = factor1_D * df2dx2;
    J.insert( 1, 1 ) = factor1 * df2dx2_D;
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
    x0 << 0.5, 1;
  }
};
