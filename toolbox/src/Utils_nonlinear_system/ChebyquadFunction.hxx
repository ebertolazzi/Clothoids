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

class ChebyquadFunction : public NonlinearSystem
{
  mutable real_type T[10];
  mutable real_type dT[10];
  real_type         integrals[10];  // Precomputed integrals of shifted Chebyshev polynomials

public:
  // Only the values N = 1, 2, 3, 4, 5, 6, 7 and 9 may be used.
  ChebyquadFunction( integer dim )
    : NonlinearSystem(
        "Chebyquad function",
        "@article {Fletcher:1965,\n"
        "  author  = {Fletcher, R.},\n"
        "  title   = {Function minimization without evaluating derivatives "
        "-- {A} review},\n"
        "  journal = {The Computer Journal},\n"
        "  year    = {1965},\n"
        "  volume  = {8},\n"
        "  pages   = {33--41},\n"
        "  doi     = {10.1093/comjnl/8.1.33},\n"
        "}\n\n"
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
        dim )
  {
    UTILS_ASSERT(
      dim > 0 && dim < 10 && dim != 8,
      "ChebyquadFunction:: dimension n = {} must be on [1,2,3,4,5,6,7,9]",
      dim );

    // Precompute integrals of shifted Chebyshev polynomials T_i^*(x) =
    // T_i(2x-1) ∫[0,1] T_i^*(x) dx = 0 for odd i, 1/(1-i^2) for even i (i > 0)
    for ( integer i = 0; i < 10; ++i )
    {
      if ( i == 0 ) { integrals[i] = 1.0; }
      else if ( i % 2 == 1 ) { integrals[i] = 0.0; }
      else
      {
        integrals[i] = 1.0 / ( 1.0 - power2( i ) );
      }
    }
  }

  void Chebyshev( real_type x ) const
  {
    // Polinomi di Chebyshev shiftati T_i^*(x) = T_i(2x-1)
    // Usiamo la ricorrenza stabile con scaling se necessario
    T[0] = 1.0;
    if ( n >= 1 )
    {
      T[1] = 2.0 * x - 1.0;

      for ( integer j = 1; j < n; ++j )
      {
        // Ricorrenza: T_{j+1}^*(x) = 2*(2x-1)*T_j^*(x) - T_{j-1}^*(x)
        T[j + 1] = 2.0 * ( 2.0 * x - 1.0 ) * T[j] - T[j - 1];

        // Controllo per evitare overflow (anche se improbabile per n <= 9)
        if ( std::abs( T[j + 1] ) > 1e100 )
        {
          // Scaling per stabilità (raro per n piccolo)
          real_type scale = 1.0 / std::abs( T[j + 1] );
          for ( integer k = 0; k <= j + 1; ++k ) { T[k] *= scale; }
        }
      }
    }
  }

  void Chebyshev_D( real_type x ) const
  {
    // Polinomi di Chebyshev shiftati e loro derivate
    T[0]  = 1.0;
    dT[0] = 0.0;

    if ( n >= 1 )
    {
      T[1]  = 2.0 * x - 1.0;
      dT[1] = 2.0;

      for ( integer j = 1; j < n; ++j )
      {
        T[j + 1]  = 2.0 * ( 2.0 * x - 1.0 ) * T[j] - T[j - 1];
        dT[j + 1] = 4.0 * T[j] + 2.0 * ( 2.0 * x - 1.0 ) * dT[j] - dT[j - 1];

        // Controllo stabilità per le derivate
        if ( std::abs( T[j + 1] ) > 1e100 || std::abs( dT[j + 1] ) > 1e100 )
        {
          real_type scale = 1.0 / std::max( std::abs( T[j + 1] ), std::abs( dT[j + 1] ) );
          for ( integer k = 0; k <= j + 1; ++k )
          {
            T[k] *= scale;
            dT[k] *= scale;
          }
        }
      }
    }
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    // Inizializza f a zero
    f.fill( 0.0 );

    // Calcola la somma dei polinomi di Chebyshev
    for ( integer j = 0; j < n; ++j )
    {
      // Verifica che x[j] sia nel dominio [0,1] per stabilità
      real_type xj = std::max( 0.0, std::min( 1.0, x( j ) ) );
      Chebyshev( xj );

      for ( integer i = 0; i < n; ++i )
      {
        // Usa T[1]...T[n] che corrispondono a T_1^*, ..., T_n^*
        f( i ) += T[i + 1];
      }
    }

    // Normalizza e sottrae l'integrale
    real_type inv_n = 1.0 / real_type( n );
    for ( integer k = 0; k < n; ++k ) { f( k ) = f( k ) * inv_n - integrals[k + 1]; }
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    real_type inv_n = 1.0 / real_type( n );
    for ( integer j = 0; j < n; ++j )
    {
      // Verifica che x[j] sia nel dominio [0,1] per stabilità
      real_type xj = std::max( 0.0, std::min( 1.0, x( j ) ) );
      Chebyshev_D( xj );

      for ( integer i = 0; i < n; ++i ) { J.insert( i, j ) = dT[i + 1] * inv_n; }
    }
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x{ x_vec[0] };
    x.resize( n );

    // Soluzioni simmetriche e ordinate in modo crescente
    switch ( n )
    {
      case 1: x << 0.5; break;
      case 2: x << 0.21132486540518711774542560974902, 0.78867513459481288225457439025098; break;
      case 3:
        x << 0.14644660940672623779957781894758, 0.50000000000000000000000000000000, 0.85355339059327376220042218105242;
        break;
      case 4:
        x << 0.10267276385411693852223453573602, 0.40620376295746005006993032654600, 0.59379623704253994993006967345400,
          0.89732723614588306147776546426398;
        break;
      case 5:
        x << 0.08375125649950906205370820817751, 0.31272929522320946720697781214099, 0.50000000000000000000000000000000,
          0.68727070477679053279302218785901, 0.91624874350049093794629179182249;
        break;
      case 6:
        x << 0.06687659094608970430820097297488, 0.28874067311944423544072680139424, 0.36668229924164763983423273308575,
          0.63331770075835236016576726691425, 0.71125932688055576455927319860576, 0.93312340905391029569179902702512;
        break;
      case 7:
        x << 0.05806914962097548214788795460948, 0.23517161235742159430747623331967, 0.33804409474004618124016345382358,
          0.50000000000000000000000000000000, 0.66195590525995381875983654617642, 0.76482838764257840569252376668033,
          0.94193085037902451785211204539052;
        break;
      case 9:
        x << 0.04420534613578276316752571608378, 0.19949067230988096428593603393306, 0.23561910847106000336990918932266,
          0.41604690789259802846598405087299, 0.50000000000000000000000000000000, 0.58395309210740197153401594912701,
          0.76438089152893999663009081067734, 0.80050932769011903571406396606694, 0.95579465386421723683247428391622;
        break;
      default:
        // Per n=8 non esiste soluzione standard - usa punti equispaziati
        for ( integer i = 0; i < n; ++i ) { x( i ) = ( i + 1.0 ) / ( n + 1.0 ); }
        break;
    }
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    for ( integer i = 0; i < n; ++i ) { x0( i ) = ( i + 1.0 ) / ( n + 1.0 ); }
  }

  virtual void bounding_box( Vector & L, Vector & U ) const override
  {
    // Restringe il bounding box a [0,1] per stabilità
    L.fill( 0.0 );
    U.fill( 1.0 );
  }
};
