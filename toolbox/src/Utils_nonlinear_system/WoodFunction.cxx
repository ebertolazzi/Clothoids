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

class WoodFunction : public NonlinearSystem
{
  using Matrix = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;

public:
  WoodFunction()
    : NonlinearSystem(
        "Wood function",
        "@book{Colville:1968,\n"
        "  author = {Colville, A.R.},\n"
        "  title  = {A Comparative Study on Nonlinear Programming Codes},\n"
        "  year   = {1968},\n"
        "  notes  = {Rep. 320-2949, New York Scientific Center}\n"
        "}\n\n"
        "@article{More:1981,\n"
        "  author  = {Mor{\'e}, Jorge J. and Garbow, Burton S. and "
        "Hillstrom, Kenneth E.},\n"
        "  title   = {Testing Unconstrained Optimization Software},\n"
        "  journal = {ACM Trans. Math. Softw.},\n"
        "  year    = {1981},\n"
        "  volume  = {7},\n"
        "  number  = {1},\n"
        "  pages   = {17--41},\n"
        "  doi     = {10.1145/355934.355936},\n"
        "}\n",
        4 )
  {
  }

  real_type t_fun( Vector const & x, integer i ) const
  {
    switch ( i )
    {
      case 0: return sqrt( 100.0 ) * ( x( 1 ) - x( 0 ) * x( 0 ) );
      case 1: return 1.0 - x( 0 );
      case 2: return sqrt( 90 ) * ( x( 3 ) - x( 2 ) * x( 2 ) );
      case 3: return 1.0 - x( 2 );
      case 4: return sqrt( 10 ) * ( x( 1 ) + x( 3 ) - 2.0 );
      case 5: return sqrt( 0.1 ) * ( x( 1 ) - x( 3 ) );
    }
    return 0;
  }

  void t_grad( Vector const & x, integer i, Vector & g ) const
  {
    g.setZero();
    switch ( i )
    {
      case 0:
        g( 0 ) = -20.00 * x( 0 );
        g( 1 ) = 10.00;
        break;
      case 1: g( 0 ) = -1; break;
      case 2:
        g( 2 ) = -6 * sqrt( 10 ) * x( 2 );
        g( 3 ) = 3 * sqrt( 10 );
        break;
      case 3: g( 2 ) = -1; break;
      case 4: g( 1 ) = g( 3 ) = sqrt( 10.0 ); break;
      case 5:
        g( 1 ) = sqrt( 0.1 );
        g( 3 ) = -sqrt( 0.1 );
        break;
    }
  }

  void t_hess( Vector const &, integer i, Matrix & h ) const
  {
    h.setZero();
    switch ( i )
    {
      case 0: h( 0, 0 ) = -20; break;
      case 1: break;
      case 2: h( 2, 2 ) = -6 * sqrt( 10.0 ); break;
      case 3: break;
      case 4: break;
      case 5: break;
    }
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    Vector g( 4 );
    f.setZero();
    for ( integer i = 0; i < 6; ++i )
    {
      real_type t = t_fun( x, i );
      t_grad( x, i, g );
      f += t * g;
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    Matrix J_full( n, n );
    J_full.setZero();
    for ( integer i = 0; i < 6; ++i )
    {
      Vector    g( 4 );
      Matrix    h( 4, 4 );
      real_type t = t_fun( x, i );
      t_grad( x, i, g );
      t_hess( x, i, h );
      J_full( 0, 0 ) += t * h( 0, 0 ) + g( 0 ) * g( 0 );
      J_full( 0, 1 ) += t * h( 0, 1 ) + g( 0 ) * g( 1 );
      J_full( 0, 2 ) += t * h( 0, 2 ) + g( 0 ) * g( 2 );
      J_full( 0, 3 ) += t * h( 0, 3 ) + g( 0 ) * g( 3 );

      J_full( 1, 1 ) += t * h( 1, 1 ) + g( 1 ) * g( 1 );
      J_full( 1, 2 ) += t * h( 1, 2 ) + g( 1 ) * g( 2 );
      J_full( 1, 3 ) += t * h( 1, 3 ) + g( 1 ) * g( 3 );

      J_full( 2, 2 ) += t * h( 2, 2 ) + g( 2 ) * g( 2 );
      J_full( 2, 3 ) += t * h( 2, 3 ) + g( 2 ) * g( 3 );

      J_full( 3, 3 ) += t * h( 3, 3 ) + g( 3 ) * g( 3 );
    }
    J_full( 1, 0 ) = J_full( 0, 1 );
    J_full( 2, 0 ) = J_full( 0, 2 );
    J_full( 3, 0 ) = J_full( 0, 3 );

    J_full( 2, 1 ) = J_full( 1, 2 );
    J_full( 3, 1 ) = J_full( 1, 3 );

    J_full( 3, 2 ) = J_full( 2, 3 );

    J.resize( n, n );
    J = J_full.sparseView();
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
    x0 << -3, -1, -3, -1;
  }
};
