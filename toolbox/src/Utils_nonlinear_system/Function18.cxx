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

class Function18 : public NonlinearSystem
{
public:
  Function18( integer neq )
    : NonlinearSystem(
        "Function 18",
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
      Vector const & x = X.segment( i, 3 );
      F( i + 0 )       = x( 0 ) * x( 1 ) - x( 2 ) * x( 2 ) - 1;
      F( i + 1 )       = x( 0 ) * x( 1 ) * x( 2 ) - x( 0 ) * x( 0 ) + x( 1 ) * x( 1 ) - 2;
      F( i + 2 )       = exp( x( 0 ) ) - exp( x( 1 ) );
    }
  }

  virtual void jacobian( Vector const & X, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer i = 0; i < n; i += 3 )
    {
      Vector const & x = X.segment( i, 3 );

      J.insert( i + 0, i + 0 ) = x( 1 );
      J.insert( i + 0, i + 1 ) = x( 0 );
      J.insert( i + 0, i + 2 ) = -2 * x( 2 );

      J.insert( i + 1, i + 0 ) = x( 1 ) * x( 2 ) - 2 * x( 0 );
      J.insert( i + 1, i + 1 ) = x( 0 ) * x( 2 ) + 2 * x( 1 );
      J.insert( i + 1, i + 2 ) = x( 0 ) * x( 1 );

      J.insert( i + 2, i + 0 ) = exp( x( 0 ) );
      J.insert( i + 2, i + 1 ) = -exp( x( 1 ) );
      J.insert( i + 2, i + 2 ) = 0;
    }
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    for ( integer i{ 0 }; i < n; i += 3 )
    {
      x0( i + 0 ) = sqrt( 2.0 );
      x0( i + 1 ) = sqrt( 2.0 );
      x0( i + 2 ) = 1;
    }
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.setZero();
  }
};
