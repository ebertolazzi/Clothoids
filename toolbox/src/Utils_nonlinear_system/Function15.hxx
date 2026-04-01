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

class Function15 : public NonlinearSystem
{
  using Matrix = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;
  typedef pair<integer, integer> INDEX;
  mutable map<INDEX, real_type>  jac_idx_vals;

public:
  Function15( integer neq )
    : NonlinearSystem(
        "Function 15",
        "@article{LaCruz:2003,\n"
        "  author    = { William {La Cruz}  and  Marcos Raydan},\n"
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
    jac_idx_vals.clear();
    jac_idx_vals[INDEX( 0, 0 )]         = 1;
    jac_idx_vals[INDEX( n - 1, n - 1 )] = 1;
    jac_idx_vals[INDEX( n - 1, n - 2 )] = 1;
    for ( integer i = 1; i < n - 1; ++i )
    {
      jac_idx_vals[INDEX( i, i - 1 )] = 1;
      jac_idx_vals[INDEX( i, i )]     = 1;
      jac_idx_vals[INDEX( i, i + 1 )] = 1;
    }
    for ( integer i = 0; i < n; ++i )
    {
      jac_idx_vals[INDEX( i, n - 5 )] = 1;
      jac_idx_vals[INDEX( i, n - 4 )] = 1;
      jac_idx_vals[INDEX( i, n - 3 )] = 1;
      jac_idx_vals[INDEX( i, n - 2 )] = 1;
      jac_idx_vals[INDEX( i, n - 1 )] = 1;
    }
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type bf = 3 * x( n - 5 ) - x( n - 4 ) - x( n - 3 ) + 0.5 * x( n - 2 ) - x( n - 1 ) + 1;
    f( 0 )       = -2 * x( 0 ) * x( 0 ) + 3 * x( 0 ) + bf;
    f( n - 1 )   = -2 * x( n - 1 ) * x( n - 1 ) + 3 * x( n - 1 ) - x( n - 2 ) + bf;
    for ( integer i = 1; i < n - 1; ++i ) f( i ) = -2 * x( i ) * x( i ) + 3 * x( i ) - x( i - 1 ) - 2 * x( i + 1 ) + bf;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    Matrix J_full( n, n );
    J_full.setZero();

    // Calcola le derivate di bf rispetto alle variabili
    // bf = 3*x(n-5) - x(n-4) - x(n-3) + 0.5 * x(n-2) - x(n-1) + 1
    // Le derivate di bf sono costanti e devono essere aggiunte a tutte le
    // equazioni

    // Prima equazione (i=0)
    J_full( 0, 0 ) = -4 * x( 0 ) + 3;
    // Aggiungi termini da bf per la prima equazione
    J_full( 0, n - 5 ) = 3;
    J_full( 0, n - 4 ) = -1;
    J_full( 0, n - 3 ) = -1;
    J_full( 0, n - 2 ) = 0.5;
    J_full( 0, n - 1 ) = -1;

    // Equazioni intermedie (i da 1 a n-2)
    for ( integer i = 1; i < n - 1; ++i )
    {
      J_full( i, i - 1 ) += -1;
      J_full( i, i ) += -4 * x( i ) + 3;
      J_full( i, i + 1 ) += -2;
      // Aggiungi termini da bf per le equazioni intermedie
      J_full( i, n - 5 ) += 3;
      J_full( i, n - 4 ) += -1;
      J_full( i, n - 3 ) += -1;
      J_full( i, n - 2 ) += 0.5;
      J_full( i, n - 1 ) += -1;
    }

    // Ultima equazione (i=n-1)
    J_full( n - 1, n - 2 ) += -1;                   // derivata di -x(n-2)
    J_full( n - 1, n - 1 ) += -4 * x( n - 1 ) + 3;  // derivata di -2*x(n-1)^2 + 3*x(n-1)
    // Aggiungi termini da bf per l'ultima equazione
    J_full( n - 1, n - 5 ) += 3;
    J_full( n - 1, n - 4 ) += -1;
    J_full( n - 1, n - 3 ) += -1;
    J_full( n - 1, n - 2 ) += 0.5;  // NOTA: += perché già inserito -1 sopra
    J_full( n - 1, n - 1 ) += -1;   // NOTA: += perché già inserito -4*x(n-1)+3 sopra

    J = J_full.sparseView();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( -1 );
  }
};
