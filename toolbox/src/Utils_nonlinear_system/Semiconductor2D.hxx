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

class Semiconductor2D : public NonlinearSystem
{
  real_type alpha;
  real_type ni;
  real_type V;
  real_type D;

public:
  Semiconductor2D()
    : NonlinearSystem(
        "2D semiconductor",
        "@techreport{Nowak:1991,\n"
        "  author = {U. Nowak and L. Weimann},\n"
        "  title  = {A Family of Newton Co des for Systems of Highly "
        "Nonlinear Equations},\n"
        "  number = {Technical Report TR-91-10 (December 1991)},\n"
        "  year   = {1991}\n"
        "}\n",
        6 )
    , alpha( 38.683 )
    , ni( 1.22E10 )
    , V( 100 )
    , D( 1E7 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = exp( alpha * ( x( 2 ) - x( 0 ) ) ) - exp( alpha * ( x( 0 ) - x( 1 ) ) ) - D / ni;
    f( 1 ) = x( 1 );
    f( 2 ) = x( 2 );
    f( 3 ) = exp( alpha * ( x( 5 ) - x( 3 ) ) ) - exp( alpha * ( x( 3 ) - x( 4 ) ) ) + D / ni;
    f( 4 ) = x( 4 ) - V;
    f( 5 ) = x( 5 ) - V;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    J.insert( 0, 0 ) = -alpha * ( exp( alpha * ( x( 0 ) - x( 1 ) ) ) + exp( alpha * ( x( 2 ) - x( 0 ) ) ) );
    J.insert( 0, 1 ) = alpha * exp( alpha * ( x( 0 ) - x( 1 ) ) );
    J.insert( 0, 2 ) = alpha * exp( alpha * ( x( 2 ) - x( 0 ) ) );

    J.insert( 1, 1 ) = 1;
    J.insert( 2, 2 ) = 1;

    J.insert( 3, 3 ) = -alpha * ( exp( alpha * ( x( 5 ) - x( 3 ) ) ) + exp( alpha * ( x( 3 ) - x( 4 ) ) ) );
    J.insert( 3, 4 ) = alpha * exp( alpha * ( x( 3 ) - x( 4 ) ) );
    J.insert( 3, 5 ) = alpha * exp( alpha * ( x( 5 ) - x( 3 ) ) );

    J.insert( 4, 4 ) = 1;
    J.insert( 5, 5 ) = 1;

    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 1 );
  }
};
