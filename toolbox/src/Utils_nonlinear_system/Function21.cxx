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

class Function21 : public NonlinearSystem
{
public:
  Function21( integer neq )
    : NonlinearSystem(
        "Function 21",
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
    check_three( n, 3 );
  }

  virtual void evaluate( Vector const & X, Vector & F ) const override
  {
    for ( integer i = 0; i < n; i += 3 )
    {
      real_type const & x = X( i + 0 );
      real_type const & y = X( i + 1 );
      real_type const & z = X( i + 2 );
      F( i + 0 )          = x * y - z * z - 1;
      F( i + 1 )          = x * y * z - x * x + y * y - 2;
      F( i + 2 )          = exp( -x ) - exp( -y );
    }
  }

  virtual void jacobian( Vector const & X, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer i = 0; i < n; i += 3 )
    {
      integer I0 = i + 0;
      integer I1 = i + 1;
      integer I2 = i + 2;

      real_type const & x = X( I0 );
      real_type const & y = X( I1 );
      real_type const & z = X( I2 );

      J.insert( I0, I0 ) = y;
      J.insert( I0, I1 ) = x;
      J.insert( I0, I2 ) = -2 * z;

      J.insert( I1, I0 ) = y * z - 2 * x;
      J.insert( I1, I1 ) = x * z + 2 * y;
      J.insert( I1, I2 ) = x * y;

      J.insert( I2, I0 ) = -exp( -x );
      J.insert( I2, I1 ) = exp( -y );
      J.insert( I2, I2 ) = 0;
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
