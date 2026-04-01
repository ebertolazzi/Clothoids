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

class Order10to11function : public NonlinearSystem
{
  real_type const SCALE;

public:
  Order10to11function()
    : NonlinearSystem(
        "Order10to11 function",
        "@article{Shacham:1972,\n"
        "  author  = {Mordechai Shacham and Ephraim Kehat},\n"
        "  title   = {An iteration method with memory for\n"
        "             the solution of a non-linear equation},\n"
        "  journal = {Chemical Engineering Science},\n"
        "  volume  = {27},\n"
        "  number  = {11},\n"
        "  pages   = {2099--2101},\n"
        "  year    = {1972},\n"
        "  doi     = {10.1016/0009-2509(72)87067-2}\n"
        "}\n",
        1 )
    , SCALE( 1e-10 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type T = x( 0 );
    f( 0 )      = SCALE * ( exp( 21000. / T ) / ( T * T ) - 1.11E11 );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type T      = x( 0 );
    real_type tmp    = 21000 / T;
    J.insert( 0, 0 ) = -SCALE * ( exp( tmp ) / power3( T ) ) * ( 2 + tmp );
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0( 0 ) = 555;
  }
};
