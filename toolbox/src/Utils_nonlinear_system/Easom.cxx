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

class Easom : public NonlinearSystem
{
public:
  Easom()
    : NonlinearSystem(
        "Easom Function",
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
    real_type arg = -power2( x( 0 ) - m_pi ) - power2( x( 1 ) - m_pi );
    f( 0 )        = cos( x( 1 ) ) * ( sin( x( 0 ) ) + 2 * cos( x( 0 ) ) * ( x( 0 ) - m_pi ) ) * exp( arg );
    f( 1 )        = cos( x( 0 ) ) * ( sin( x( 1 ) ) + 2 * cos( x( 1 ) ) * ( x( 1 ) - m_pi ) ) * exp( arg );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    real_type arg     = -power2( x( 0 ) - m_pi ) - power2( x( 1 ) - m_pi );
    real_type dargdx1 = -2 * ( x( 0 ) - m_pi );
    real_type dargdx2 = -2 * ( x( 1 ) - m_pi );

    real_type factor = cos( x( 1 ) ) * ( sin( x( 0 ) ) - cos( x( 0 ) ) * dargdx1 );
    real_type dfdx1  = cos( x( 1 ) ) * ( cos( x( 0 ) ) + sin( x( 0 ) ) * dargdx1 + 2 * cos( x( 0 ) ) );
    real_type dfdx2  = -sin( x( 1 ) ) * ( sin( x( 0 ) ) - cos( x( 0 ) ) * dargdx1 );

    J.insert( 0, 0 ) = ( dfdx1 + factor * dargdx1 ) * exp( arg );
    J.insert( 0, 1 ) = ( dfdx2 + factor * dargdx2 ) * exp( arg );

    factor = cos( x( 0 ) ) * ( sin( x( 1 ) ) - cos( x( 1 ) ) * dargdx2 );
    dfdx1  = -sin( x( 0 ) ) * ( sin( x( 1 ) ) - cos( x( 1 ) ) * dargdx2 );
    dfdx2  = cos( x( 0 ) ) * ( cos( x( 1 ) ) + sin( x( 1 ) ) * dargdx2 + 2 * cos( x( 1 ) ) );

    J.insert( 1, 0 ) = ( dfdx1 + factor * dargdx1 ) * exp( arg );
    J.insert( 1, 1 ) = ( dfdx2 + factor * dargdx2 ) * exp( arg );

    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( m_pi );
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 0.5, 1;
  }
};
