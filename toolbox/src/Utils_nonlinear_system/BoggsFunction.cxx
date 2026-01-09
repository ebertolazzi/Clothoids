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

class BoggsFunction : public NonlinearSystem
{
public:
  BoggsFunction()
    : NonlinearSystem(
        "Boggs function",
        "@article{Boggs:1971,\n"
        "  author = {Boggs, P.},\n"
        "  title  = {The Solution of Nonlinear Systems of Equations\n"
        "            by A-Stable Integration Techniques},\n"
        "  journal = {SIAM Journal on Numerical Analysis},\n"
        "  volume  = {8},\n"
        "  number  = {4},\n"
        "  pages   = {767--785},\n"
        "  year    = {1971},\n"
        "  doi     = {10.1137/0708071},\n"
        "}\n",
        2 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = power2( x( 0 ) ) - x( 1 ) + 1;
    f( 1 ) = x( 0 ) - cos( m_pi_2 * x( 1 ) );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    J.insert( 0, 0 ) = 2 * x( 0 );
    J.insert( 0, 1 ) = -1;
    J.insert( 1, 0 ) = 1;
    J.insert( 1, 1 ) = m_pi_2 * sin( m_pi_2 * x( 1 ) );
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 0, 1;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 1, 0;
  }
};
