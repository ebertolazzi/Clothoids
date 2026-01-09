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

class McKinnon : public NonlinearSystem
{
  real_type const tau, theta, phi;

public:
  McKinnon()
    : NonlinearSystem(
        "McKinnon function (stable)",
        "@article{McKinnon:1998,\n"
        "  author  = {McKinnon, K.},\n"
        "  title   = {Convergence of the Nelder--Mead Simplex\n"
        "             Method to a Nonstationary Point},\n"
        "  journal = {SIAM Journal on Optimization},\n"
        "  volume  = {9},\n"
        "  number  = {1},\n"
        "  pages   = {148-158},\n"
        "  year    = {1998},\n"
        "  doi     = {10.1137/S1052623496303482},\n"
        "}\n",
        2 )
    , tau( 2.0 )
    , theta( 6.0 )
    , phi( 60.0 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type const eps    = 1e-12;  // tolleranza per x0 ~ 0
    real_type       x0     = x( 0 );
    real_type       x0_abs = std::abs( x0 );

    if ( x0_abs < eps )
    {
      // serie di Taylor attorno a 0: f0 ~ theta * tau * x0^(tau-1)
      f( 0 ) = ( x0 >= 0 ? 1.0 : phi ) * theta * tau * std::pow( eps, tau - 1 ) * ( x0 / eps );
    }
    else
    {
      f( 0 ) = ( x0 >= 0 ? 1.0 : phi ) * theta * tau * std::pow( x0_abs, tau - 1 );
      if ( x0 < 0 ) f( 0 ) = -f( 0 );  // segno corretto per x0<0
    }

    f( 1 ) = 1.0 + 2.0 * x( 1 );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    real_type const eps    = 1e-12;
    real_type       x0     = x( 0 );
    real_type       x0_abs = std::abs( x0 );

    if ( x0_abs < eps )
    {
      // derivata tramite limite / Taylor
      J.insert( 0, 0 ) = ( x0 >= 0 ? 1.0 : phi ) * theta * tau * ( tau - 1 ) * std::pow( eps, tau - 2 );
    }
    else
    {
      real_type sign   = ( x0 >= 0 ? 1.0 : -1.0 );
      J.insert( 0, 0 ) = sign * ( x0 >= 0 ? 1.0 : phi ) * theta * tau * ( tau - 1 ) * std::pow( x0_abs, tau - 2 );
    }

    J.insert( 1, 1 ) = 2.0;
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0_ref = x_vec[0];
    x0_ref.resize( n );
    x0_ref << 0, -0.5;  // soluzione nota
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0_ref = x_vec[0];
    x0_ref.resize( n );
    x0_ref.fill( 1.0 );
  }
};
