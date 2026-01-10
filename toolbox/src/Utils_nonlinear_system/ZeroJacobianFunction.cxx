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

class ZeroJacobianFunction : public NonlinearSystem
{
public:
  ZeroJacobianFunction( integer neq )
    : NonlinearSystem(
        "Zero Jacobian Function (same as function 27)",
        "@article{LaCruz:2003,\n"
        "  author    = {William {La Cruz}  and  Marcos Raydan},\n"
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
    check_min_equations( n, 1 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = x.dot( x );
    for ( integer i = 1; i < n; ++i ) f( i ) = -2 * x( 0 ) * x( i );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer i = 0; i < n; ++i ) J.insert( 0, i ) = 2 * x( i );
    for ( integer i = 1; i < n; ++i )
    {
      J.insert( i, 0 ) = -2 * x( i );
      J.insert( i, i ) = -2 * x( 0 );
    }
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    real_type bf = ( 1.0 / 60.0 - 100.0 / ( 6.0 * n ) ) * ( 1.0 / 60.0 - 50.0 / 6.0 );
    x0.fill( bf );
    x0( 0 ) = 100.0 * ( n - 100.0 ) / n;
  }
};
