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

class JennrichAndSampsonFunction : public NonlinearSystem
{
public:
  JennrichAndSampsonFunction()
    : NonlinearSystem(
        "Jennrich and Sampson function",
        "@article{Jennrich:1968,\n"
        "  author  = { Jennrich, R. I. and Sampson, P. F.},\n"
        "  title   = {Application of Stepwise Regression to Non-Linear "
        "Estimation},\n"
        "  journal = {Technometrics},\n"
        "  year    = {1968},\n"
        "  volume  = {10},\n"
        "  number  = {1},\n"
        "  pages   = {63--72},\n"
        "  doi     = {10.1080/00401706.1968.10490535},\n"
        "}\n\n"
        "@article{More:1981,\n"
        "  author  = {Mor{\'e}, Jorge J. and Garbow, Burton S. and "
        "Hillstrom, Kenneth E.},\n"
        "  title   = {Testing Unconstrained Optimization Software},\n"
        "  journal = {ACM Trans. Math. Softw.},\n"
        "  year    = {1981},\n"
        "  volume  = {7},\n"
        "  number  = {1},\n"
        "  pages   = {17--41},\n"
        "  doi     = {10.1145/355934.355936},\n"
        "}\n",
        2 )
  {
  }

  virtual void evaluate( Vector const & x_in, Vector & f ) const override
  {
    real_type x   = x_in[0];
    real_type y   = x_in[1];
    real_type t1  = exp( x );
    real_type t2  = exp( y );
    real_type t3  = 4.0 - t1 - t2;
    real_type t7  = exp( 2.0 * x );
    real_type t9  = exp( 2.0 * y );
    real_type t10 = 6.0 - t7 - t9;
    f( 0 )        = -2.0 * t3 * t1 - 4.0 * t10 * t7;
    f( 1 )        = -2.0 * t3 * t2 - 4.0 * t10 * t9;
  }

  virtual void jacobian( Vector const & x_in, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type x      = x_in[0];
    real_type y      = x_in[1];
    real_type t1     = exp( x );
    real_type t2     = t1 * t1;
    real_type t4     = exp( y );
    real_type t5     = 4.0 - t1 - t4;
    real_type t9     = exp( 2.0 * x );
    real_type t10    = t9 * t9;
    real_type t13    = exp( 2.0 * y );
    real_type t14    = 6.0 - t9 - t13;
    real_type t22    = 2.0 * t4 * t1 + 8.0 * t13 * t9;
    real_type t23    = t4 * t4;
    real_type t27    = t13 * t13;
    J.insert( 0, 0 ) = 2.0 * t2 - 2.0 * t5 * t1 + 8.0 * t10 - 8.0 * t14 * t9;
    J.insert( 0, 1 ) = t22;
    J.insert( 1, 0 ) = t22;
    J.insert( 1, 1 ) = 2.0 * t23 - 2.0 * t5 * t4 + 8.0 * t27 - 8.0 * t14 * t13;
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 0.5609475731693775630903801021735698602305366752976206558112028362911293020590355239361946223920 );
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 0.3, 0.4;
  }
};
