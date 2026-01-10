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

class DiscreteIntegralEquationFunction : public NonlinearSystem
{
public:
  DiscreteIntegralEquationFunction( integer neq )
    : NonlinearSystem(
        "Discrete integral equation function",
        "@article{More:1981,\n"
        "  author  = {Mor{\'e}, Jorge J. and Garbow, Burton S. and "
        "Hillstrom, Kenneth E.},\n"
        "  title   = {Testing Unconstrained Optimization Software},\n"
        "  journal = {ACM Trans. Math. Softw.},\n"
        "  year    = {1981},\n"
        "  volume  = {7},\n"
        "  number  = {1},\n"
        "  pages   = {17--41},\n"
        "  doi     = {10.1145/355934.355936},\n"
        "}\n",
        neq )
  {
    check_min_equations( n, 2 );
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type h = 1 / real_type( n + 1 );

    for ( integer k{ 0 }; k < n; ++k )
    {
      real_type tk   = ( k + 1.0 ) / ( n + 1.0 );
      real_type sum1 = 0;
      for ( integer j{ 0 }; j < k; ++j )
      {
        real_type tj = ( j + 1 ) * h;
        sum1 += tj * power3( x( j ) + tj + 1 );
      }
      real_type sum2 = 0;
      for ( integer j{ k }; j < n; ++j )
      {
        real_type tj = ( j + 1 ) * h;
        sum2 += ( 1 - tj ) * power3( x( j ) + tj + 1 );
      }
      f( k ) = x( k ) + ( h / 2 ) * ( ( 1 - tk ) * sum1 + tk * sum2 );
    }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    for ( integer k = 0; k < n; ++k )
    {
      real_type tk = real_type( k + 1 ) / real_type( n + 1 );
      for ( integer j = 0; j < n; ++j )
      {
        real_type tj    = real_type( j + 1 ) / real_type( n + 1 );
        real_type temp1 = power2( x( j ) + tj + 1 );
        real_type temp2 = min( tk, tj ) - tj * tk;
        real_type val   = 1.5 * temp2 * temp1 / real_type( n + 1 );
        if ( j == k ) val += 1;
        J.insert( k, j ) = val;
      }
    }
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.clear();
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    switch ( n )
    {
      case 2: x0 << -0.1282467630337316243876390987485189444347, -0.1592675672446408630494624980879210125394; break;
      case 3:
        x0 << -0.1048174251268790134252362144215285375609, -0.1627022933732565189436998826409546303353,
          -0.1458503921843392374798322333967325583654;
        break;
      case 4:
        x0 << -0.08764857853130447170049222985553575439651, -0.1477703373755026679125035808546396389734,
          -0.1686201922723071170660805480770337174359, -0.1308164496352389909774613554127782652712;
        break;
      case 5:
        x0 << -0.07502212929232048088309048985370741564914, -0.1319762103521906431069193729290020520734,
          -0.1648487719093373112876960905846226308130, -0.1646646802158007321519060418145486053884,
          -0.1174176516841936157849999649609128856169;
        break;
      case 10:
        x0 << -0.0431649825187648705768864628080, -0.0754165336858920839547720879022,
          -0.0815771565353868815338983704309, -0.114485714380529287243408349831, -0.140973576862596679632277275305,
          -0.159908696181983122325572374408, -0.169877202312774918976188661992, -0.169089983781208351844136654063,
          -0.155249535221831821946902537001, -0.125355891678934989400407791540;
        break;
      default: x_vec.clear(); break;
    }
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    for ( integer k{ 0 }; k < n; ++k ) x0( k ) = real_type( ( k + 1 ) * ( k - n ) ) / power2( n + 1.0 );
  }
};
