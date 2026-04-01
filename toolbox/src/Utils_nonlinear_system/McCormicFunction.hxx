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

class McCormicFunction : public NonlinearSystem
{
public:
  McCormicFunction()
    : NonlinearSystem(
        "McCormic function",
        "MVF - Multivariate Test Functions Library in C for Unconstrained "
        "Global Optimization\n"
        "Ernesto P. Adorio Department of Mathematics U.P. Diliman\n"
        "ernesto.adorio@gmail.com eadorio@yahoo.com, 2005\n",
        2 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type t2 = cos( x( 0 ) + x( 1 ) );
    real_type t3 = 2.0 * x( 0 );
    real_type t4 = 2.0 * x( 1 );
    f( 0 )       = t2 + t3 - t4 - 3.0 / 2.0;
    f( 1 )       = t2 - t3 + t4 + 5.0 / 2.0;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type t2     = sin( x( 0 ) + x( 1 ) );
    real_type t3     = -t2 + 2.0;
    real_type t4     = -t2 - 2.0;
    J.insert( 0, 0 ) = t3;
    J.insert( 0, 1 ) = t4;
    J.insert( 1, 0 ) = t4;
    J.insert( 1, 1 ) = t3;
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << m_pi / 3 + 0.5, m_pi / 3 - 0.5;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 2 );
    auto & x0{ x_vec[0] };
    auto & x1{ x_vec[1] };
    x0.resize( n );
    x1.resize( n );
    x0 << -1.5, 4;
    x1 << -3, 4;
  }
};
