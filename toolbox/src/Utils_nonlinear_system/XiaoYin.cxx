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

class XiaoYin1 : public NonlinearSystem
{
public:
  XiaoYin1()
    : NonlinearSystem(
        "XiaoYin example 3",
        "@article{XiaoYin:2015,\n"
        "  author    = {Xiaoyong Xiao and Hongwei Yin},\n"
        "  title     = {A new class of methods with higher order of "
        "convergence for solving systems of nonlinear equations},\n"
        "  Journal   = {Applied Mathematics and Computation},\n"
        "  Number    = {264},\n"
        "  Pages     = {300–-309},\n"
        "  Publisher = {Elsevier},\n"
        "  Year      = {2015},\n"
        "  Doi       = {10.1016/j.amc.2015.04.094},\n"
        "}\n",
        8 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    for ( integer i = 0; i < n; ++i ) f( i ) = x( i ) * log( 1 + x( ( i + 1 ) % n ) ) - 1;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer i = 0; i < n; ++i )
    {
      J.insert( i, i )             = log( 1 + x( ( i + 1 ) % n ) );
      J.insert( i, ( i + 1 ) % n ) = x( i ) / ( 1 + x( ( i + 1 ) % n ) );
    }
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 1.5 );
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class XiaoYin2 : public NonlinearSystem
{
public:
  XiaoYin2()
    : NonlinearSystem(
        "XiaoYin example 4",
        "@article{XiaoYin:2015,\n"
        "  author    = {Xiaoyong Xiao and Hongwei Yin},\n"
        "  title     = {A new class of methods with higher order of "
        "convergence for solving systems of nonlinear equations},\n"
        "  Journal   = {Applied Mathematics and Computation},\n"
        "  Number    = {264},\n"
        "  Pages     = {300–-309},\n"
        "  Publisher = {Elsevier},\n"
        "  Year      = {2015},\n"
        "  Doi       = {10.1016/j.amc.2015.04.094},\n"
        "}\n",
        16 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    for ( integer i = 0; i < n; ++i ) f( i ) = x( i ) * sin( x( ( i + 1 ) % n ) ) - 1;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer i = 0; i < n; ++i )
    {
      J.insert( i, i )             = sin( x( ( i + 1 ) % n ) );
      J.insert( i, ( i + 1 ) % n ) = x( i ) * cos( x( ( i + 1 ) % n ) );
    }
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( -0.85 );
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class XiaoYin3 : public NonlinearSystem
{
public:
  XiaoYin3()
    : NonlinearSystem(
        "XiaoYin example 5",
        "@article{XiaoYin:2015,\n"
        "  author    = {Xiaoyong Xiao and Hongwei Yin},\n"
        "  title     = {A new class of methods with higher order of "
        "convergence for solving systems of nonlinear equations},\n"
        "  Journal   = {Applied Mathematics and Computation},\n"
        "  Number    = {264},\n"
        "  Pages     = {300–-309},\n"
        "  Publisher = {Elsevier},\n"
        "  Year      = {2015},\n"
        "  Doi       = {10.1016/j.amc.2015.04.094},\n"
        "}\n",
        35 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    for ( integer i = 0; i < n; ++i )
    {
      real_type x1  = x( i );
      real_type xi1 = x( ( i + 1 ) % n );
      f( i )        = x1 * xi1 - exp( -x1 ) - exp( -xi1 );
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer i = 0; i < n; ++i )
    {
      real_type x1                 = x( i );
      real_type xi1                = x( ( i + 1 ) % n );
      J.insert( i, i )             = xi1 + exp( -x1 );
      J.insert( i, ( i + 1 ) % n ) = x1 + exp( -xi1 );
    }
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 1.2 );
  }
};
