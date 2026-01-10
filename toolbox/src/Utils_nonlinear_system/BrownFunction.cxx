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

class BrownFunction : public NonlinearSystem
{
public:
  BrownFunction()
    : NonlinearSystem(
        "Brown function",
        "@article{Qi:2006,\n"
        "  author  = {Qi, H. and Sun, D.},\n"
        "  title   = {A Quadratically Convergent Newton Method for\n"
        "             Computing the Nearest Correlation Matrix},\n"
        "  journal = {SIAM Journal on Matrix Analysis and Applications},\n"
        "  volume  = {28},\n"
        "  number  = {2},\n"
        "  pages   = {360--385},\n"
        "  year    = {2006},\n"
        "  doi     = {10.1137/050624509},\n"
        "}\n",
        2 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = power2( x( 0 ) ) - ( x( 1 ) + 1 );
    f( 1 ) = power2( x( 0 ) - 2 ) + power2( x( 1 ) - 0.5 ) - 1;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    J.insert( 0, 0 ) = 2 * x( 0 );
    J.insert( 0, 1 ) = -1;
    J.insert( 1, 0 ) = 2 * x( 0 ) - 4;
    J.insert( 1, 1 ) = 2 * x( 1 ) - 1;

    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 1.067346085806689713408597312807010489088, 0.1392276668868614404836249880514108793336;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 0.1, 2;
  }
};
