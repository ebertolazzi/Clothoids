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

class YixunShi1 : public NonlinearSystem
{
public:
  YixunShi1()
    : NonlinearSystem(
        "Shi, Yixun Problem N.3",
        "@article{YixunShi,\n"
        "  Author  = {Shi, Yixun},\n"
        "  Title   = {A globalization procedure for solving "
        "nonlinear systems of equations},\n"
        "  Journal = {Numerical Algorithms},\n"
        "  Number  = {2},\n"
        "  Pages   = {273--286},\n"
        "  Volume  = {12},\n"
        "  Year    = {1996},\n"
        "  Doi     = {10.1007/BF02142807},\n"
        "}\n",
        3 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = 3 * x( 0 ) - cos( x( 1 ) * x( 2 ) ) - 0.5;
    f( 1 ) = x( 0 ) * x( 0 ) - 625 * x( 1 ) * x( 1 );
    f( 2 ) = exp( -x( 0 ) * x( 1 ) ) + 20 * x( 2 ) + ( 10 * m_pi - 3 ) / 3;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type x0     = x( 0 );
    real_type x1     = x( 1 );
    real_type x2     = x( 2 );
    real_type t2     = sin( x1 * x2 );
    real_type t8     = exp( -x0 * x1 );
    J.insert( 0, 0 ) = 3;
    J.insert( 0, 1 ) = x2 * t2;
    J.insert( 0, 2 ) = x1 * t2;
    J.insert( 1, 0 ) = 2 * x0;
    J.insert( 1, 1 ) = -1250 * x1;
    J.insert( 2, 0 ) = -x1 * t8;
    J.insert( 2, 1 ) = -x0 * t8;
    J.insert( 2, 2 ) = 20;
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0( 0 ) = 0.4999816893677071161061094979125234746691;
    x0( 1 ) = -0.1999926757470828464424437991650093898676e-1;
    x0( 2 ) = -0.5241012469638837054562573412192983462072;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 3 );
    auto & x0{ x_vec[0] };
    auto & x1{ x_vec[1] };
    auto & x2{ x_vec[2] };
    x0.resize( n );
    x1.resize( n );
    x2.resize( n );
    x0.fill( 1 );
    x1 << 0, 1e-6, 0;
    x1.setZero();
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class YixunShi2 : public NonlinearSystem
{
public:
  YixunShi2()
    : NonlinearSystem(
        "Shi, Yixun Problem N.4",
        "@article{YixunShi,\n"
        "  Author  = {Shi, Yixun},\n"
        "  Title   = {A globalization procedure for solving "
        "nonlinear systems of equations},\n"
        "  Journal = {Numerical Algorithms},\n"
        "  Number  = {2},\n"
        "  Pages   = {273--286},\n"
        "  Volume  = {12},\n"
        "  Year    = {1996},\n"
        "  Doi     = {10.1007/BF02142807},\n"
        "}\n",
        100 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f = x.array() * ( 3 - 2 * x.array() ) + 1 + x( 49 ) / 2;
    f( 0 ) -= 2 * x( 1 );
    f( 99 ) -= x( 98 );
    for ( integer idx = 1; idx < 99; ++idx ) f( idx ) -= 2 * x( idx + 1 ) + x( idx - 1 );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    J.insert( 0, 0 )  = 3 - 4 * x( 0 );
    J.insert( 0, 1 )  = -2;
    J.insert( 0, 49 ) = 0.5;
    for ( integer idx = 1; idx < 99; ++idx )
    {
      real_type tmp1 = -1;
      real_type tmp2 = 3 - 4 * x( idx );
      real_type tmp3 = -2;
      if ( idx == 48 )
      {
        tmp1 += 0.5;  // (48,49)
      }
      else if ( idx == 49 )
      {
        tmp2 += 0.5;  // (49,49)
      }
      else if ( idx == 50 )
      {
        tmp3 += 0.5;  // (50,49)
      }
      else
      {
        J.insert( idx, 49 ) = 0.5;
      }
      J.insert( idx, idx - 1 ) = tmp1;
      J.insert( idx, idx )     = tmp2;
      J.insert( idx, idx + 1 ) = tmp3;
    }
    J.insert( 99, 49 ) = 0.5;
    J.insert( 99, 98 ) = -1;
    J.insert( 99, 99 ) = 3 - 4 * x( 99 );
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 4 );
    auto & x0{ x_vec[0] };
    auto & x1{ x_vec[1] };
    auto & x2{ x_vec[2] };
    auto & x3{ x_vec[3] };
    x0.resize( n );
    x1.resize( n );
    x2.resize( n );
    x3.resize( n );
    x0.fill( -1 );
    x1.fill( 0.04 );
    x2.fill( 0.044285 );
    x3.fill( 1 );
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class YixunShi3 : public NonlinearSystem
{
public:
  YixunShi3()
    : NonlinearSystem(
        "Shi, Yixun Problem N.5",
        "@article{YixunShi,\n"
        "  Author  = {Shi, Yixun},\n"
        "  Title   = {A globalization procedure for solving "
        "nonlinear systems of equations},\n"
        "  Journal = {Numerical Algorithms},\n"
        "  Number  = {2},\n"
        "  Pages   = {273--286},\n"
        "  Volume  = {12},\n"
        "  Year    = {1996},\n"
        "  Doi     = {10.1007/BF02142807},\n"
        "}\n",
        100 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f = 1 - x( 99 ) + 0.5 * x( 98 ) - x( 97 ) - x( 96 ) + 3 * x( 95 ) + x.array() * ( 3 - 2 * x.array() );
    f( 0 ) -= 2 * x( 1 );
    f( 99 ) -= x( 98 );
    for ( integer idx = 1; idx < 99; ++idx ) f( idx ) -= 2 * x( idx + 1 ) + x( idx - 1 );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer idx = 0; idx < 100; ++idx )
    {
      real_type tmp0 = 3;    // 95
      real_type tmp1 = -1;   // 96
      real_type tmp2 = -1;   // 97
      real_type tmp3 = 0.5;  // 98
      real_type tmp4 = -1;   // 99
      if ( idx == 0 )
      {
        J.insert( 0, 0 ) = 3 - 4 * x( 0 );
        J.insert( 0, 1 ) = -2;
      }
      else if ( idx < 94 )
      {
        J.insert( idx, idx - 1 ) = -1;
        J.insert( idx, idx )     = 3 - 4 * x( idx );
        J.insert( idx, idx + 1 ) = -2;
      }
      else if ( idx == 94 )
      {
        tmp0 -= 2;  // (94,95)
        J.insert( 94, 93 ) = -1;
        J.insert( 94, 94 ) = 3 - 4 * x( idx );
      }
      else if ( idx == 95 )
      {
        tmp0 += 3 - 4 * x( idx );  // (95,95)
        tmp1 -= 2;                 // (95,96)
        J.insert( 95, 94 ) = -1;
      }
      else if ( idx == 96 )
      {
        tmp0 -= 1;
        tmp1 += 3 - 4 * x( idx );
        tmp2 -= 2;
      }
      else if ( idx == 97 )
      {
        tmp1 -= 1;
        tmp2 += 3 - 4 * x( idx );
        tmp3 -= 2;
      }
      else if ( idx == 98 )
      {
        tmp2 -= 1;
        tmp3 += 3 - 4 * x( idx );
        tmp4 -= 2;
      }
      else if ( idx == 99 )
      {
        tmp3 -= 1;
        tmp4 += 3 - 4 * x( idx );
      }
      J.insert( idx, 95 ) = tmp0;  // 95
      J.insert( idx, 96 ) = tmp1;  // 96
      J.insert( idx, 97 ) = tmp2;  // 97
      J.insert( idx, 98 ) = tmp3;  // 98
      J.insert( idx, 99 ) = tmp4;  // 99
    }
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 4 );
    auto & x0{ x_vec[0] };
    auto & x1{ x_vec[1] };
    auto & x2{ x_vec[2] };
    auto & x3{ x_vec[3] };
    x0.resize( n );
    x1.resize( n );
    x2.resize( n );
    x3.resize( n );
    x0.fill( -1 );
    x1.fill( 0.04 );
    x2.fill( 0.043283 );
    x3.fill( 1 );
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class YixunShi4 : public NonlinearSystem
{
  real_type h;

public:
  YixunShi4()
    : NonlinearSystem(
        "Shi, Yixun Problem N.6 (Singular Broyden)",
        "@article{YixunShi,\n"
        "  Author  = {Shi, Yixun},\n"
        "  Title   = {A globalization procedure for solving "
        "nonlinear systems of equations},\n"
        "  Journal = {Numerical Algorithms},\n"
        "  Number  = {2},\n"
        "  Pages   = {273--286},\n"
        "  Volume  = {12},\n"
        "  Year    = {1996},\n"
        "  Doi     = {10.1007/BF02142807},\n"
        "}\n",
        100 )
    , h( 0.5 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f = ( 3 - h * x.array() ) * x.array() + 1;
    f( 0 ) -= 2 * x( 1 );
    f( 99 ) -= x( 98 );
    for ( integer idx = 1; idx < 99; ++idx ) f( idx ) -= 2 * x( idx + 1 ) + x( idx - 1 );
    f = f.array() * f.array();
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type tmp    = ( 3 - h * x( 0 ) ) * x( 0 ) - 2 * x( 1 ) + 1;
    real_type tmp_1  = 3 - 2 * h * x( 0 );
    J.insert( 0, 0 ) = 2 * tmp * tmp_1;
    J.insert( 0, 1 ) = -4 * tmp;
    for ( integer idx = 1; idx < 99; ++idx )
    {
      tmp                      = ( 3 - h * x( idx ) ) * x( idx ) - x( idx - 1 ) - 2 * x( idx + 1 ) + 1;
      tmp_1                    = 3 - 2 * h * x( idx );
      J.insert( idx, idx - 1 ) = -2 * tmp;
      J.insert( idx, idx )     = 2 * tmp * tmp_1;
      J.insert( idx, idx + 1 ) = -4 * tmp;
    }
    tmp                = ( 3 - h * x( 99 ) ) * x( 99 ) - x( 98 ) + 1;
    tmp_1              = 3 - 2 * h * x( 99 );
    J.insert( 99, 98 ) = -2 * tmp;
    J.insert( 99, 99 ) = 2 * tmp * tmp_1;
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 6 );
    auto & x0{ x_vec[0] };
    auto & x1{ x_vec[1] };
    auto & x2{ x_vec[2] };
    auto & x3{ x_vec[3] };
    auto & x4{ x_vec[4] };
    auto & x5{ x_vec[5] };
    x0.resize( n );
    x1.resize( n );
    x2.resize( n );
    x3.resize( n );
    x4.resize( n );
    x5.resize( n );

    x0.fill( -1 );
    x1.fill( -0.001 );
    x2.setZero();
    x3.fill( 0.1 );
    x4.fill( 0.17 );
    x5.fill( 1 );
  }
};
