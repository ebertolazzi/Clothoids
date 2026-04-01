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

class HanSunHan : public NonlinearSystem
{
  real_type rr[99];

public:
  HanSunHan()
    : NonlinearSystem(
        "Han-Sun-Han-SAMPAJO 2005 function test",
        "@article{Han:2005,\n"
        "  author  = {Qiaoming Han and Wenyu Sun and Jiye Han "
        "and Raimudo J. B. Sampaio},\n"
        "  title   = {An adaptive conic trust-region method "
        "for unconstrained optimization},\n"
        "  journal = {Optimization Methods and Software},\n"
        "  year    = {2005},\n"
        "  volume  = {20},\n"
        "  number  = {6},\n"
        "  pages   = {665--677},\n"
        "  doi = {10.1080/10556780410001697677}\n"
        "}\n",
        2 )
  {
    for ( integer i = 0; i < 99; ++i )
    {
      real_type arg = ( i + 1 ) / 100.0;
      rr[i]         = pow( -50.0 * log( arg ), 2.0 / 3.0 ) + 25.0;
    }
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type t1  = x( 0 ) * x( 0 );
    real_type t3  = t1 * t1;
    real_type t7  = x( 1 ) * x( 1 );
    real_type t12 = exp( t7 );
    f( 0 )        = ( 8 * t3 * t1 * x( 0 ) ) + ( 2 * x( 0 ) ) + 2 * x( 0 ) * t7;
    f( 1 )        = ( 2 * t1 * x( 1 ) ) + 2 * x( 1 ) * t12;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type t1     = x( 0 ) * x( 0 );
    real_type t2     = t1 * t1;
    real_type t5     = x( 1 ) * x( 1 );
    real_type t9     = 4 * x( 0 ) * x( 1 );
    real_type t11    = exp( t5 );
    J.insert( 0, 0 ) = ( 56 * t2 * t1 ) + 2 + 2 * t5;
    J.insert( 0, 1 ) = J.insert( 1, 0 ) = t9;
    J.insert( 1, 1 )                    = ( 2 * t1 ) + 2 * t11 + 4 * t5 * t11;
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 5, 3;
  }
};
