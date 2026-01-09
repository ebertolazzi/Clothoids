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

class HanbookFunction : public NonlinearSystem
{
  mutable real_type sum1;
  mutable real_type sum2;

public:
  HanbookFunction( integer neq )
    : NonlinearSystem(
        "Hanbook Function",
        "@techreport{Raydan:2004,\n"
        "  author = {William La Cruz and Jose Mario Martinez and Marcos "
        "Raydan},\n"
        "  title  = {Spectral residual method without gradient\n"
        "             information for solving large-scale nonlinear\n"
        "             systems of equations: Theory and experiments},\n"
        "  number = {Technical Report RT-04-08},\n"
        "  year   = {2004}\n"
        "}\n\n"
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

  void sum( Vector const & x ) const
  {
    sum1 = sum2 = 0;
    for ( integer i = 0; i < n; ++i )
    {
      real_type xm = x( i ) - 1;
      sum1 += xm;
      sum2 += xm * xm;
    }
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    sum( x );
    real_type S12 = sin( sum1 + sum2 );
    real_type S1  = 2 * sin( sum1 );
    for ( integer i = 0; i < n; ++i ) f( i ) = 0.05 * ( x( i ) - 1 ) + ( 2 + 4 * ( x( i ) - 1 ) ) * S12 + S1;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    sum( x );
    real_type S12 = sin( sum1 + sum2 );
    real_type C12 = cos( sum1 + sum2 );
    real_type C1  = 2 * cos( sum1 );
    for ( integer i = 0; i < n; ++i )
    {
      for ( integer j = 0; j < n; ++j )
      {
        real_type tmp = ( 2 + 4 * ( x( i ) - 1 ) ) * ( 2 * x( j ) - 1 ) * C12 + C1;
        if ( i == j ) tmp += 0.05 + 4 * S12;
        J.insert( i, j ) = tmp;
      }
    }
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 5 );
  }
};
