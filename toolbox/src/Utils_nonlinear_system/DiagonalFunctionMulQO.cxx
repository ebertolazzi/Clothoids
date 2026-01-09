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

class DiagonalFunctionMulQO : public NonlinearSystem
{
public:
  DiagonalFunctionMulQO( integer neq )
    : NonlinearSystem(
        "Diagonal Functions Multiplied by quasi-orthogonal matrix",
        "@article{Gasparo:2000,\n"
        "  Author    = {Maria Grazia Gasparo},\n"
        "  Title     = {A nonmonotone hybrid method for nonlinear "
        "systems},\n"
        "  Journal   = {Optimization Methods and Software},\n"
        "  Number    = {2},\n"
        "  Pages     = {79--94},\n"
        "  Publisher = {Taylor & Francis},\n"
        "  Volume    = {13},\n"
        "  Year      = {2000},\n"
        "  Doi       = {10.1080/10556780008805776},\n"
        "}\n",
        neq )
  {
    check_three( n, 3 );
  }

  virtual void evaluate( Vector const & X, Vector & F ) const override
  {
    for ( integer i = 0; i < n; i += 3 )
    {
      real_type x0 = X( i + 0 );
      real_type x1 = X( i + 1 );
      real_type x2 = X( i + 2 );
      F( i + 0 )   = x0 * ( 0.6 + 1.6 * x0 * x0 ) + x1 * ( 9.6 - 7.2 * x1 ) - 4.8;
      F( i + 1 )   = 0.48 * x0 + x1 * ( -4.32 + x1 * ( 3.24 - 0.72 * x1 ) ) + x2 * ( 0.2 * x2 * x2 - 1 ) + 2.16;
      F( i + 2 )   = x2 * ( 1.25 - 0.25 * x2 * x2 );
    }
  }

  virtual void jacobian( Vector const & X, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer i = 0; i < n; i += 3 )
    {
      real_type x0 = X( i + 0 );
      real_type x1 = X( i + 1 );
      real_type x2 = X( i + 2 );

      J.insert( i + 0, i + 0 ) = 0.6 + 4.8 * x0 * x0;
      J.insert( i + 0, i + 1 ) = 9.6 - 14.4 * x1;

      J.insert( i + 1, i + 0 ) = 0.48;
      J.insert( i + 1, i + 1 ) = -4.32 + x1 * ( 6.48 - 2.16 * x1 );
      J.insert( i + 1, i + 2 ) = 0.6 * x2 * x2 - 1;

      J.insert( i + 2, i + 2 ) = 1.25 - 0.75 * x2 * x2;
    }
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    for ( integer i = 0; i < n; i += 3 )
    {
      x0( i + 0 ) = -1;
      x0( i + 1 ) = 0.5;
      x0( i + 2 ) = -1;
    }
  }
};
