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

class BraninRCOS : public NonlinearSystem
{
  real_type const a, d, e, b, c, ff;

public:
  BraninRCOS()
    : NonlinearSystem(
        "BraninRCOS",
        "@book{brent2013,\n"
        "  author    = {Brent, R.P.},\n"
        "  title     = {Algorithms for Minimization Without Derivatives},\n"
        "  isbn      = {9780486143682},\n"
        "  series    = {Dover Books on Mathematics},\n"
        "  year      = {2013},\n"
        "  publisher = {Dover Publications}\n"
        "}\n",
        2 )
    , a( 1.0 )
    , d( 6.0 )
    , e( 10.0 )
    , b( 5.1 / ( 4.0 * m_pi * m_pi ) )
    , c( 5.0 / m_pi )
    , ff( 1.0 / ( 8.0 * m_pi ) )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type x1 = x( 0 );
    real_type x2 = x( 1 );
    f( 0 )       = 2.0 * a * ( x2 - b * x1 * x1 + c * x1 - d ) * ( c - 2 * b * x1 ) - e * ( 1 - ff ) * sin( x1 );
    f( 1 )       = 2.0 * a * ( x2 - b * x1 * x1 + c * x1 - d );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    real_type x1 = x( 0 );
    real_type x2 = x( 1 );

    J.insert( 0, 0 ) = 2 * a * power2( c - 2 * b * x1 ) - 4 * a * b * ( x2 - b * x1 * x1 + c * x1 - d ) -
                       e * ( 1 - ff ) * cos( x1 );
    J.insert( 0, 1 ) = J.insert( 1, 0 ) = 2 * a * ( c - 2 * b * x1 );
    J.insert( 1, 1 )                    = 2 * a;

    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
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

    x0 << -m_pi, 12.275;
    x1 << m_pi, 2.275;
    x2 << 9.4247779607693797153879301498385086525915081981254, 2.475;
    x3 << 0, 6;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << -1, 1;
  }
};
