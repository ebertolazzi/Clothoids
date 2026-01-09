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

class BrownAndConteFunction : public NonlinearSystem
{
  real_type const cst;

public:
  BrownAndConteFunction()
    : NonlinearSystem(
        "Brown and Conte function",
        "@inproceedings{Brown:1967,\n"
        "  author    = {Brown, Kenneth M. and Conte, Samuel D.},\n"
        "  title     = {The Solution of Simultaneous Nonlinear "
        "Equations},\n"
        "  booktitle = {Proceedings of the 1967 22Nd National "
        "Conference},\n"
        "  series    = {ACM '67},\n"
        "  year      = {1967},\n"
        "  pages     = {111--114},\n"
        "  doi       = {10.1145/800196.805981},\n"
        "  acmid     = {805981},\n"
        "  publisher = {ACM},\n"
        "}\n",
        2 )
    , cst( 1 - 1 / ( 4 * m_pi ) )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = ( sin( x( 0 ) * x( 1 ) ) - x( 1 ) / ( 2 * m_pi ) - x( 0 ) ) / 2;
    f( 1 ) = cst * ( exp( 2 * x( 0 ) ) - m_e ) + x( 1 ) * m_e * m_1_pi - 2 * m_e * x( 0 );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    J.insert( 0, 0 ) = ( x( 1 ) * cos( x( 0 ) * x( 1 ) ) - 1 ) / 2;
    J.insert( 0, 1 ) = x( 0 ) * cos( x( 0 ) * x( 1 ) ) / 2 - 0.25 / m_pi;
    J.insert( 1, 0 ) = 2 * ( cst * exp( 2 * x( 0 ) ) - m_e );
    J.insert( 1, 1 ) = m_e * m_1_pi;
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 0.5, m_pi;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 0.6, 3;
  }
};
