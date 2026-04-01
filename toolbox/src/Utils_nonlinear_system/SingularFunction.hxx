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

class SingularFunction : public NonlinearSystem
{
public:
  SingularFunction( integer neq )
    : NonlinearSystem(
        "Singular Function",
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
    check_min_equations( n, 2 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 )     = power3( x( 0 ) ) / 3 + power2( x( 1 ) ) / 2;
    f( n - 1 ) = power2( x( n - 1 ) ) * ( ( n / 3.0 ) * x( n - 1 ) - 0.5 );
    for ( integer i = 1; i < n - 1; ++i )
      f( i ) = power2( x( i ) ) * ( ( ( i + 1 ) / 3.0 ) * x( i ) - 0.5 ) + 0.5 * power2( x( i + 1 ) );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    J.insert( 0, 0 )         = power2( x( 0 ) );
    J.insert( 0, 1 )         = x( 1 );
    J.insert( n - 1, n - 1 ) = ( n * x( n - 1 ) - 1 ) * x( n - 1 );
    for ( integer i = 1; i < n - 1; ++i )
    {
      J.insert( i, i )     = ( ( i + 1 ) * x( i ) - 1 ) * x( i );
      J.insert( i, i + 1 ) = x( i + 1 );
    }
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
