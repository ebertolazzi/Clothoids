/*\
 |
 |  Author:
 |    Enrico Bertolazzi
 |    University of Trento
 |    Department of Industrial Engineering
 |    Via Sommarive 9, I-38123, Povo, Trento, Italy
 |    email: enrico.bertolazzi@unitn.it
\*/

class ArtificialTestOfNowakAndWeimann : public NonlinearSystem
{
public:
  ArtificialTestOfNowakAndWeimann()
    : NonlinearSystem(
        "Artificial Test of Nowak and Weimann",
        "@techreport{Nowak1991,\n"
        "  author = {U. Nowak and L. Weimann},\n"
        "  title  = {A Family of Newton Codes for Systems of "
        "Highly Nonlinear Equations},\n"
        "  number = {Technical Report TR-91-10},\n"
        "  year   = {1991}\n"
        "}\n",
        2 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = exp( power2( x( 0 ) ) + power2( x( 1 ) ) ) - 3;
    f( 1 ) = x( 0 ) + x( 1 ) - sin( 3 * ( x( 0 ) + x( 1 ) ) );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type tmp1   = 2 * exp( power2( x( 0 ) ) + power2( x( 1 ) ) );
    J.insert( 0, 0 ) = x( 0 ) * tmp1;
    J.insert( 0, 1 ) = x( 1 ) * tmp1;

    real_type tmp2   = 1 - 3 * cos( 3 * ( x( 0 ) + x( 1 ) ) );
    J.insert( 1, 0 ) = tmp2;
    J.insert( 1, 1 ) = tmp2;
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 0.81, 0.82;
  }
};
;
