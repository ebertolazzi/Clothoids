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

class BadlyScaledAugmentedPowellFunction : public NonlinearSystem
{
  real_type phi( real_type t ) const
  {
    if ( t <= -1 )
      return t / 2 - 2;
    else if ( t >= 2 )
      return t / 2 + 2;
    else
      return ( -1924 + t * ( 4551 + t * ( 888 - t * 592 ) ) ) / 1998;
  }

  real_type phi_1( real_type t ) const
  {
    if ( t <= -1 )
      return 0.5;
    else if ( t >= 2 )
      return 0.5;
    else
      return ( 4551 + t * ( 2 * 888 - t * 3 * 592 ) ) / 1998;
  }

public:
  BadlyScaledAugmentedPowellFunction( integer neq )
    : NonlinearSystem(
        "Badly scaled augmented Powell’s function",
        "@article{Gasparo:2000,\n"
        "  Author    = {Maria Grazia Gasparo},\n"
        "  Title     = {A nonmonotone hybrid method for "
        "nonlinear systems},\n"
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
    for ( integer k = 0; k < n; k += 3 )
    {
      Vector const & x = X.segment( k, 3 );
      F( k + 0 )       = 10000 * ( x( 0 ) * x( 1 ) ) - 1.0;
      F( k + 1 )       = exp( -x( 1 ) ) + exp( -x( 0 ) ) - 1.0001;
      F( k + 2 )       = phi( x( 2 ) );
    }
  }

  virtual void jacobian( Vector const & X, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer k = 0; k < n; k += 3 )
    {
      Vector const & x = X.segment( k, 3 );

      // F(k+0) = 10000 * (x(0) * x(1)) - 1.0
      J.insert( k + 0, k + 0 ) = 10000 * x( 1 );  // derivata rispetto a x(0)
      J.insert( k + 0, k + 1 ) = 10000 * x( 0 );  // derivata rispetto a x(1)

      // F(k+1) = exp(-x(1)) + exp(-x(0)) - 1.0001
      J.insert( k + 1, k + 0 ) = -exp( -x( 0 ) );  // derivata di exp(-x(0)) rispetto a x(0)
      J.insert( k + 1, k + 1 ) = -exp( -x( 1 ) );  // derivata di exp(-x(1)) rispetto a x(1)

      // F(k+2) = phi(x(2))
      J.insert( k + 2, k + 2 ) = phi_1( x( 2 ) );  // derivata di phi(x(2)) rispetto a x(2)
    }
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    for ( integer k{ 0 }; k < n; k += 3 )
    {
      x0( k + 0 ) = 0.109815932969981745568376164563E-4;
      x0( k + 1 ) = 9.10614673986652401094671049032;
      x0( k + 2 ) = 0.3998810580736440979319618294548679254646;
    }
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 2 );
    auto & x0{ x_vec[0] };
    auto & x1{ x_vec[1] };
    x0.resize( n );
    x1.resize( n );
    for ( integer k{ 0 }; k < n; k += 3 )
    {
      x0( k + 0 ) = 0;
      x0( k + 1 ) = 1;
      x0( k + 2 ) = -4;
    }
    for ( integer k = 0; k < n; k += 3 )
    {
      x1( k + 0 ) = 1e-3;
      x1( k + 1 ) = 18;
      x1( k + 2 ) = 1;
    }
  }
};
