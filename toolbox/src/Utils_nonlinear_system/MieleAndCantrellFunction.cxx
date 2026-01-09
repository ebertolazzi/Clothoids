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

class MieleAndCantrellFunction : public NonlinearSystem
{
  using Matrix = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;

public:
  MieleAndCantrellFunction()
    : NonlinearSystem(
        "Miele and Cantrell function",
        "@article{Grippo:1991,\n"
        "  author  = {Grippo, L. and Lampariello, F. and Lucidi, S.},\n"
        "  title   = {A Class of Nonmonotone Stabilization Methods\n"
        "             in Unconstrained Optimization},\n"
        "  journal = {Numer. Math.},\n"
        "  year    = {1991},\n"
        "  volume  = {59},\n"
        "  number  = {1},\n"
        "  pages   = {779--805},\n"
        "  doi     = {10.1007/BF01385810},\n"
        "}\n",
        4 )
  {
  }

private:
  static constexpr real_type eps = 1e-12;  // soglia per x0 vicino a zero

  inline real_type safe_power( real_type x, int n ) const
  {
    if ( std::abs( x ) < eps ) return 0.0;
    return std::pow( x, n );
  }

  /* Gradient contributions */
  inline void add_grad1( const Vector & x, Vector & g ) const
  {
    const real_type e  = std::exp( x( 0 ) );
    const real_type d  = e - x( 1 );
    const real_type d3 = d * d * d;
    const real_type c  = 4 * d3;
    g( 0 ) += c * e;
    g( 1 ) -= c;
  }

  inline void add_grad2( const Vector & x, Vector & g ) const
  {
    const real_type d = x( 1 ) - x( 2 );
    g( 1 ) += 600 * d * d * d * d * d;
    g( 2 ) -= 600 * d * d * d * d * d;
  }

  inline void add_grad3( const Vector & x, Vector & g ) const
  {
    const real_type z = x( 2 ) - x( 3 );
    const real_type t = std::tan( z );
    g( 2 ) += 4 * t * t * t * ( 1 + t * t );
    g( 3 ) -= 4 * t * t * t * ( 1 + t * t );
  }

  /* Hessian contributions */
  inline void add_hess1( const Vector & x, Matrix & h ) const
  {
    const real_type e   = std::exp( x( 0 ) );
    const real_type t   = e - x( 1 );
    const real_type t2  = t * t;
    const real_type t2e = t2 * e;
    h( 0, 0 ) += 4 * t2e * ( 4 * e - 1 );
    h( 0, 1 ) -= 12 * t2e;
    h( 1, 0 ) -= 12 * t2e;
    h( 1, 1 ) += 12 * t2;
  }

  inline void add_hess2( const Vector & x, Matrix & h ) const
  {
    const real_type d   = x( 1 ) - x( 2 );
    const real_type d4  = d * d * d * d;
    const real_type tmp = 3000 * d4;
    h( 1, 1 ) += tmp;
    h( 1, 2 ) -= tmp;
    h( 2, 1 ) -= tmp;
    h( 2, 2 ) += tmp;
  }

  inline void add_hess3( const Vector & x, Matrix & h ) const
  {
    const real_type z    = x( 2 ) - x( 3 );
    const real_type t    = std::tan( z );
    const real_type t2   = t * t;
    const real_type sec2 = 1 + t2;
    const real_type tmp  = 4 * t2 * sec2 * ( 3 + 5 * t2 + 2 * t2 * t2 );
    h( 2, 2 ) += tmp;
    h( 2, 3 ) -= tmp;
    h( 3, 2 ) -= tmp;
    h( 3, 3 ) += tmp;
  }

public:
  void evaluate( const Vector & x, Vector & f ) const override
  {
    f.setZero();
    f( 0 ) = 56 * safe_power( x( 0 ), 6 );  // protezione numerica per x0 ~ 0
    add_grad2( x, f );
    add_grad3( x, f );
    add_grad1( x, f );
  }

  void jacobian( const Vector & x, SparseMatrix & J ) const override
  {
    Matrix h( 4, 4 );
    h.setZero();
    h( 0, 0 ) = 336 * safe_power( x( 0 ), 5 );  // d²/dx0² con protezione
    add_hess2( x, h );
    add_hess3( x, h );
    add_hess1( x, h );

    J.resize( n, n );
    J.setZero();
    for ( integer i = 0; i < 4; i++ )
      for ( integer j = 0; j < 4; j++ )
        if ( h( i, j ) != 0.0 ) J.insert( i, j ) = h( i, j );
    J.makeCompressed();
  }

  void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    x_vec[0].resize( 4 );
    x_vec[0] << 0, 1, 1, 1;
  }

  void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 2 );
    x_vec[0].resize( 4 );
    x_vec[1].resize( 4 );
    x_vec[0] << 10, -10, -10, -10;
    x_vec[1] << 1, 2, 2, 2;
  }
};
