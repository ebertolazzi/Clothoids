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

class ShenYpma5 : public NonlinearSystem
{
public:
  ShenYpma5()
    : NonlinearSystem(
        "Shen-Ypma Example N.5",
        "@article{,\n"
        "  Author = {Yun-Qiu Shen and Tjalling J. Ypma},\n"
        "  Doi = {10.1016/j.apnum.2004.09.029},\n"
        "  Journal = {Applied Numerical Mathematics},\n"
        "  Number = {2},\n"
        "  Pages = {256 - 265},\n"
        "  Title = {Newton's method for singular nonlinear equations using "
        "approximate left and right nullspaces of the Jacobian},\n"
        "  Volume = {54},\n"
        "  Year = {2005},\n"
        "}\n",
        2 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type x0 = x( 0 );
    real_type x1 = x( 1 );
    f( 0 )       = x0 * x0 * ( 1 - x0 * x1 ) + x1 * x1;
    f( 1 )       = x0 * x0 + x1 * x1 * ( 3 * x0 - 2 );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type x0     = x( 0 );
    real_type x1     = x( 1 );
    J.insert( 0, 0 ) = x0 * ( 2 - 3 * x0 * x1 );
    J.insert( 0, 1 ) = 2 * x1 - x0 * x0 * x0;
    J.insert( 1, 0 ) = 2 * x0 + 3 * x1 * x1;
    J.insert( 1, 1 ) = x1 * ( 6 * x0 - 4 );
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
    x_vec.resize( 2 );
    auto & x0{ x_vec[0] };
    auto & x1{ x_vec[1] };
    x0.resize( n );
    x1.resize( n );

    x0.fill( 0.02 );
    x1 << 6.2989, 3.7048;
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class ShenYpma7 : public NonlinearSystem
{
public:
  ShenYpma7()
    : NonlinearSystem(
        "Shen-Ypma Example N.7",
        "@article{,\n"
        "  Author = {Yun-Qiu Shen and Tjalling J. Ypma},\n"
        "  Doi = {10.1016/j.apnum.2004.09.029},\n"
        "  Journal = {Applied Numerical Mathematics},\n"
        "  Number = {2},\n"
        "  Pages = {256 - 265},\n"
        "  Title = {Newton's method for singular nonlinear equations using "
        "approximate left and right nullspaces of the Jacobian},\n"
        "  Volume = {54},\n"
        "  Year = {2005},\n"
        "}\n",
        2 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type x0 = x( 0 );
    real_type x1 = x( 1 );
    f( 0 )       = x0 * x0 - x1 * x1;
    f( 1 )       = 3 * ( x0 * x0 - x1 * x1 );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type x0     = x( 0 );
    real_type x1     = x( 1 );
    J.insert( 0, 0 ) = 2 * x0;
    J.insert( 0, 1 ) = -2 * x1;
    J.insert( 1, 0 ) = 6 * x0;
    J.insert( 1, 1 ) = -6 * x1;
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
    x_vec.resize( 2 );
    auto & x0{ x_vec[0] };
    auto & x1{ x_vec[1] };
    x0.resize( n );
    x1.resize( n );

    x0 << 0.05, 0.04;
    x1 << 1, 2;
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class ShenYpma8 : public NonlinearSystem
{
public:
  ShenYpma8()
    : NonlinearSystem(
        "Shen-Ypma Example N.8",
        "@article{,\n"
        "  Author = {Yun-Qiu Shen and Tjalling J. Ypma},\n"
        "  Doi = {10.1016/j.apnum.2004.09.029},\n"
        "  Journal = {Applied Numerical Mathematics},\n"
        "  Number = {2},\n"
        "  Pages = {256 - 265},\n"
        "  Title = {Newton's method for singular nonlinear equations using "
        "approximate left and right nullspaces of the Jacobian},\n"
        "  Volume = {54},\n"
        "  Year = {2005},\n"
        "}\n",
        5 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type x1 = x( 0 );
    real_type x2 = x( 1 );
    real_type x3 = x( 2 );
    real_type x4 = x( 3 );
    real_type x5 = x( 4 );
    f( 0 )       = x1 + x2 + x3 * x3 + x4 * x4 + x5 * x5 - 2;
    f( 1 )       = x1 - x2 + x3 * x3 + x4 * x4 + x5 * x5;
    f( 2 )       = -x3 * x3 + x4 * x4 + x5 * x5;
    f( 3 )       = x3 * x3 - x4 * x4 + x5 * x5;
    f( 4 )       = x3 * x3 + x4 * x4 - x5 * x5;
    ;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    // real_type x1 = x(0);
    // real_type x2 = x(1);
    real_type x3 = x( 2 );
    real_type x4 = x( 3 );
    real_type x5 = x( 4 );

    J.resize( n, n );
    J.setZero();

    J.insert( 0, 0 ) = 1;
    J.insert( 0, 1 ) = 1;
    J.insert( 0, 2 ) = 2 * x3;
    J.insert( 0, 3 ) = 2 * x4;
    J.insert( 0, 4 ) = 2 * x5;

    J.insert( 1, 0 ) = 1;
    J.insert( 1, 1 ) = -1;
    J.insert( 1, 2 ) = 2 * x3;
    J.insert( 1, 3 ) = 2 * x4;
    J.insert( 1, 4 ) = 2 * x5;

    J.insert( 2, 0 ) = 0;
    J.insert( 2, 1 ) = 0;
    J.insert( 2, 2 ) = -2 * x3;
    J.insert( 2, 3 ) = 2 * x4;
    J.insert( 2, 4 ) = 2 * x5;

    J.insert( 3, 0 ) = 0;
    J.insert( 3, 1 ) = 0;
    J.insert( 3, 2 ) = 2 * x3;
    J.insert( 3, 3 ) = -2 * x4;
    J.insert( 3, 4 ) = 2 * x5;

    J.insert( 4, 0 ) = 0;
    J.insert( 4, 1 ) = 0;
    J.insert( 4, 2 ) = 2 * x3;
    J.insert( 4, 3 ) = 2 * x4;
    J.insert( 4, 4 ) = -2 * x5;

    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 1 );
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 1.02, 1.02, 0.02, 0.02, 0.02;
  }
};
