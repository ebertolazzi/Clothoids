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

class ComplementaryFunction : public NonlinearSystem
{
public:
  ComplementaryFunction( integer neq )
    : NonlinearSystem(
        "Complementary Function",
        "@article{LaCruz:2006,\n"
        "  title   = {Spectral Residual Method without Gradient "
        "Information\n"
        "             for Solving Large-Scale Nonlinear Systems of "
        "Equations},\n"
        "  author  = {William La Cruz and Jos\\`e Mario Mart\\`\\inez and "
        "Marcos Raydan},\n"
        "  journal = {Mathematics of Computation},\n"
        "  year    = {2006},\n"
        "  volume  = {75},\n"
        "  number  = {255},\n"
        "  pages   = {1429--1448},\n"
        "  publisher = {American Mathematical Society},\n"
        "}\n\n"
        "@article{LaCruz:2003,\n"
        "  author    = { William {La Cruz}  and  Marcos Raydan},\n"
        "  title     = {Nonmonotone Spectral Methods for Large-Scale "
        "Nonlinear Systems},\n"
        "  journal   = {Optimization Methods and Software},\n"
        "  year      = {2003},\n"
        "  volume    = {18},\n"
        "  number    = {5},\n"
        "  pages     = {583--599},\n"
        "  publisher = {Taylor & Francis},\n"
        "  doi       = {10.1080/10556780310001610493},\n"
        "}\n",
        neq )
  {
    check_even( n, 2 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    for ( integer i = 0; i < n; i += 2 )
    {
      real_type t1 = x( i ) * x( i );
      real_type t2 = exp( x( i ) );
      real_type t4 = 1.0 / n;
      real_type t6 = power2( t2 * x( i ) - t4 );
      real_type t8 = sqrt( t1 + t6 );
      f( i )       = t8 - ( t2 + 1.0 ) * x( i ) + t4;
    }
    for ( integer i = 1; i < n; i += 2 )
    {
      real_type t1 = x( i ) * x( i );
      real_type t3 = sin( x( i ) );
      real_type t4 = exp( x( i ) );
      real_type t6 = power2( 3.0 * x( i ) + t3 + t4 );
      real_type t8 = sqrt( t1 + t6 );
      f( i )       = t8 - 4 * x( i ) - t3 - t4;
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer i = 0; i < n; i += 2 )
    {
      real_type t1     = x( i ) * x( i );
      real_type t2     = exp( x( i ) );
      real_type t3     = t2 * x( i );
      real_type t5     = t3 - 1.0 / n;
      real_type t6     = t5 * t5;
      real_type t8     = sqrt( t1 + t6 );
      J.insert( i, i ) = ( x( i ) + ( t2 + t3 ) * t5 ) / t8 - 1.0 - t2 - t3;
    }
    for ( integer i = 1; i < n; i += 2 )
    {
      real_type t1     = x( i ) * x( i );
      real_type t3     = sin( x( i ) );
      real_type t4     = exp( x( i ) );
      real_type t5     = 3.0 * x( i ) + t3 + t4;
      real_type t6     = t5 * t5;
      real_type t8     = sqrt( t1 + t6 );
      real_type t10    = cos( x( i ) );
      J.insert( i, i ) = ( x( i ) + ( 3.0 + t10 + t4 ) * t5 ) / t8 - 4.0 - t10 - t4;
    }
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 0.5 );
  }
};
