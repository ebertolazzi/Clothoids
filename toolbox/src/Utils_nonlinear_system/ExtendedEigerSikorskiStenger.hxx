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

class ExtendedEigerSikorskiStenger : public NonlinearSystem
{
public:
  ExtendedEigerSikorskiStenger()
    : NonlinearSystem(
        "Extended Eiger-Sikorski-Stenger Function",
        "@article{Eiger:1984,\n"
        "  author  = {Eiger, A. and Sikorski, K. and Stenger, F.},\n"
        "  title   = {A Bisection Method for Systems of Nonlinear "
        "Equations},\n"
        "  journal = {ACM Trans. Math. Softw.},\n"
        "  year    = {1984},\n"
        "  volume  = {10},\n"
        "  number  = {4},\n"
        "  pages   = {367--377},\n"
        "  doi     = {10.1145/2701.2705},\n"
        "}\n\n"
        "@article{Kearfott:1987,\n"
        "  author  = {Kearfott, R. Baker},\n"
        "  title   = {Some Tests of Generalized Bisection},\n"
        "  journal = {ACM Trans. Math. Softw.},\n"
        "  year    = {1987},\n"
        "  volume  = {13},\n"
        "  number  = {3},\n"
        "  pages   = {197--220},\n"
        "  doi     = {10.1145/29380.29862},\n"
        "}\n",
        9 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    for ( integer i = 0; i < n - 1; ++i ) f( i ) = power2( x( i ) - 0.1 ) + x( i + 1 ) - 0.1;
    f( n - 1 ) = power2( x( n - 1 ) - 0.1 ) + x( 0 ) - 0.1;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer i = 0; i < n - 1; ++i )
    {
      J.insert( i, i )     = 2 * ( x( i ) - 0.1 );
      J.insert( i, i + 1 ) = 1;
    }
    J.insert( n - 1, n - 1 ) = 2 * ( x( n - 1 ) - 0.1 );
    J.insert( n - 1, 0 )     = 1;
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( -2000 );
  }
};
