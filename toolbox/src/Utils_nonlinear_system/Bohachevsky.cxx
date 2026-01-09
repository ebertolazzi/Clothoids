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

#define Bohachevsky_BIBTEX                                                     \
  "@book{Michalewicz:1996,\n"                                                  \
  "  author = {Michalewicz, Zbigniew},\n"                                      \
  "  title = {Genetic Algorithms + Data Structures = Evolution Programs (3rd " \
  "Ed.)},\n"                                                                   \
  "  year = {1996},\n"                                                         \
  "  isbn = {3-540-60676-9},\n"                                                \
  "  publisher = {Springer-Verlag},\n"                                         \
  "  address = {Berlin, Heidelberg},\n"                                        \
  "}\n\n"                                                                      \
  "@book{brent2002algorithms,\n"                                               \
  "  author={Brent, R.P.},\n"                                                  \
  "  title={Algorithms for Minimization Without Derivatives},\n"               \
  "  year={2002},\n"                                                           \
  "  isbn={9780486419985},\n"                                                  \
  "  series={Dover Books on Mathematics},\n"                                   \
  "  publisher={Dover Publications}\n"                                         \
  "}\n"

class BohachevskyN1 : public NonlinearSystem
{
public:
  BohachevskyN1() : NonlinearSystem( "BohachevskyN1", Bohachevsky_BIBTEX, 2 ) {}

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type x1 = x( 0 );
    real_type x2 = x( 1 );
    f( 0 )       = 2.0 * x1 + 0.9 * m_pi * sin( 3.0 * m_pi * x1 );
    f( 1 )       = 4.0 * x2 + 1.6 * m_pi * sin( 4.0 * m_pi * x2 );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type x1 = x( 0 );
    real_type x2 = x( 1 );

    J.insert( 0, 0 ) = 2.0 + 2.7 * ( m_pi * m_pi ) * cos( 3.0 * m_pi * x1 );
    J.insert( 0, 1 ) = 0.0;

    J.insert( 1, 0 ) = 0.0;
    J.insert( 1, 1 ) = 4.0 + 6.4 * ( m_pi * m_pi ) * cos( 4.0 * m_pi * x2 );
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

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class BohachevskyN2 : public NonlinearSystem
{
public:
  BohachevskyN2() : NonlinearSystem( "BohachevskyN2", Bohachevsky_BIBTEX, 2 ) {}

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type x1 = x( 0 );
    real_type x2 = x( 1 );
    f( 0 )       = 2 * x1 + 0.9 * m_pi * sin( 3 * m_pi * x1 ) * cos( 4 * m_pi * x2 );
    f( 1 )       = 4 * x2 + 1.2 * m_pi * cos( 3 * m_pi * x1 ) * sin( 4 * m_pi * x2 );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type x1     = x( 0 );
    real_type x2     = x( 1 );
    real_type pi2    = m_pi * m_pi;
    real_type A      = pi2 * sin( 3 * m_pi * x1 ) * sin( 4 * m_pi * x2 );
    real_type B      = pi2 * cos( 3 * m_pi * x1 ) * cos( 4 * m_pi * x2 );
    J.insert( 0, 0 ) = 2 + 2.7 * B;
    J.insert( 0, 1 ) = -3.6 * A;
    J.insert( 1, 0 ) = -3.6 * A;
    J.insert( 1, 1 ) = 4 + 4.8 * B;
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
    x0 << 0.6, 1.3;
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class BohachevskyN3 : public NonlinearSystem
{
public:
  BohachevskyN3() : NonlinearSystem( "BohachevskyN3", Bohachevsky_BIBTEX, 2 ) {}

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type x1 = x( 0 );
    real_type x2 = x( 1 );
    f( 0 )       = 2.0 * x1 + 0.9 * m_pi * sin( 3 * m_pi * x1 );
    f( 1 )       = 4.0 * x2 - 4.0 * m_pi * sin( 4 * m_pi * x2 );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type x1     = x( 0 );
    real_type x2     = x( 1 );
    real_type pi2    = m_pi * m_pi;
    J.insert( 0, 0 ) = 2.0 + 2.7 * pi2 * cos( 3 * m_pi * x1 );
    J.insert( 1, 1 ) = 4.0 - 16.0 * pi2 * cos( 4 * m_pi * x2 );
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
    x0 << 0.5, 1.0;
  }
};
