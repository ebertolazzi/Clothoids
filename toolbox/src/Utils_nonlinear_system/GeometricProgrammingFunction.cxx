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

class GeometricProgrammingFunction : public NonlinearSystem
{
  using Matrix = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;

public:
  GeometricProgrammingFunction( integer neq )
    : NonlinearSystem(
        "Geometric Programming Function",
        "@techreport{Raydan:2004,\n"
        "  author = {William La Cruz and Jose Mario Martinez and Marcos "
        "Raydan},\n"
        "  title  = {Spectral residual method without gradient\n"
        "             information for solving large-scale nonlinear\n"
        "             systems of equations: Theory and experiments},\n"
        "  number = {Technical Report RT-04-08},\n"
        "  year   = {2004}\n"
        "}\n",
        neq )
  {
    check_min_equations( n, 2 );
  }

  // ============================================================
  // f(x)
  // ============================================================
  void evaluate( Vector const & x, Vector & f ) const override
  {
    f.fill( -1.0 );

    // precompute log(x)
    Vector logx = x.array().log();

    for ( integer t = 1; t <= 5; ++t )
    {
      real_type t1 = 0.2 * t;
      real_type t2 = t1 - 1;

      for ( integer i = 0; i < n; ++i )
      {
        real_type acc = t1;

        for ( integer k = 0; k < n; ++k )
        {
          if ( k != i )
            acc *= exp( t1 * logx( k ) );
          else
            acc *= exp( t2 * logx( k ) );
        }

        f( i ) += acc;
      }
    }

    // final product term
    for ( integer i = 0; i < n; ++i )
    {
      real_type prod = 1.0;
      for ( integer k = 0; k < n; ++k )
      {
        if ( k != i ) prod *= x( k );
      }
      f( i ) += prod;
    }
  }

  // ============================================================
  // J(x)
  // ============================================================
  void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    Matrix Jfull( n, n );
    Jfull.setZero();

    Vector logx = x.array().log();

    for ( integer t = 1; t <= 5; ++t )
    {
      real_type t1 = 0.2 * t;
      real_type t2 = t1 - 1;
      real_type t3 = t2 - 1;

      for ( integer i = 0; i < n; ++i )
      {
        for ( integer j = 0; j < n; ++j )
        {
          real_type acc = 1.0;

          if ( i == j )
          {
            // diagonal: ∂f_i/∂x_i
            acc *= t1 * t2 * exp( t3 * logx( i ) );

            for ( integer k = 0; k < n; ++k )
            {
              if ( k != i ) acc *= exp( t1 * logx( k ) );
            }
          }
          else
          {
            // off diagonal: ∂f_i/∂x_j
            for ( integer k = 0; k < n; ++k )
            {
              if ( k == i || k == j )
                acc *= t1 * std::exp( t2 * logx( k ) );
              else
                acc *= std::exp( t1 * logx( k ) );
            }
          }
          Jfull( i, j ) += acc;
        }
      }
    }
    // ============================================================
    // Derivative of product term Π_{k≠i} x_k
    // ============================================================

    for ( integer i = 0; i < n; ++i )
    {
      for ( integer j = 0; j < n; ++j )
      {
        if ( i == j ) continue;

        real_type prod = 1.0;

        // derivative wrt x_j eliminates x_j from product
        for ( integer k = 0; k < n; ++k )
        {
          if ( k == i || k == j ) continue;
          prod *= x( k );
        }

        Jfull( i, j ) += prod;
      }
    }

    // sparse output
    J = Jfull.sparseView();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    // soluzione solo per n > 2
    switch ( n )
    {
      case 3: x0.fill( 0.1998358397137673483536787376895800631111 ); break;
      case 4: x0.fill( 0.4504094438079780141262271051664676993245 ); break;
      case 5: x0.fill( 0.5752676827564249696803041659800640074052 ); break;
      case 6: x0.fill( 0.6535160063654996023801965308489436831018 ); break;
      case 7: x0.fill( 0.7073521857908597934305906573214506518189 ); break;
      case 8: x0.fill( 0.7466921578931958109831452261293280849340 ); break;
      case 9: x0.fill( 0.7767040434991580564723357591126807540423 ); break;
      case 10: x0.fill( 0.8003562243848881877566691646546685244726 ); break;
      case 50: x0.fill( 0.9618792339620107524542660085692669885314 ); break;
      case 100: x0.fill( 0.9810472935851392047929164385524051070425 ); break;
      default: x_vec.clear();
    }
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 1 );
  }

  virtual void check_if_admissible( Vector const & x ) const override
  {
    for ( integer i = 0; i < n - 1; ++i ) UTILS_ASSERT( x( i ) > 0, "x[{}] = {} must be > 0", i, x( i ) );
  }
};
