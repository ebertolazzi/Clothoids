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

class BurdenAndFaires : public NonlinearSystem
{
public:
  BurdenAndFaires()
    : NonlinearSystem(
        "Burden and Faires example 1",
        "@book{burden2005,\n"
        "  author    = {Burden, R. and Faires, J.},\n"
        "  title     = {Numerical Analysis},\n"
        "  year      = {2005},\n"
        "  pages     = {597--640},\n"
        "  publisher = {Thomson Brooks/Cole}\n"
        "}\n",
        3 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = 3 * x( 0 ) - cos( x( 1 ) * x( 2 ) ) - 0.5;
    f( 1 ) = x( 0 ) * x( 0 ) - 81 * power2( x( 1 ) + 0.1 ) + sin( x( 2 ) ) + 1.06;
    f( 2 ) = exp( -x( 0 ) * x( 1 ) ) + 20 * x( 2 ) + ( 10 * m_pi - 3 ) / 3;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    J.insert( 0, 0 ) = 3;
    J.insert( 0, 1 ) = sin( x( 1 ) * x( 2 ) ) * x( 2 );
    J.insert( 0, 2 ) = sin( x( 1 ) * x( 2 ) ) * x( 1 );

    J.insert( 1, 0 ) = 2 * x( 0 );
    J.insert( 1, 1 ) = -162 * ( x( 1 ) + 0.1 );
    J.insert( 1, 2 ) = cos( x( 2 ) );

    J.insert( 2, 0 ) = -exp( -x( 0 ) * x( 1 ) ) * x( 1 );
    J.insert( 2, 1 ) = -exp( -x( 0 ) * x( 1 ) ) * x( 0 );
    J.insert( 2, 2 ) = 20;

    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 0.1, 0.1, -0.1;
  }
};
