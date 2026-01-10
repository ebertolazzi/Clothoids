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

class CompressibilityFactorFromTheRKequation : public NonlinearSystem
{
  real_type Q, r;

public:
  CompressibilityFactorFromTheRKequation()
    : NonlinearSystem(
        "Compressibility factor from the RK equation",
        "@book{Cutlip:2007,\n"
        "  author    = {Cutlip, Michael and Shacham, Mordechai},\n"
        "  title     = {Problem Solving in Chemical and Biochemical "
        "Engineering\n"
        "               with Polymath,\\texttrademark Excel,\n"
        "               and Matlab\\textregistered, Second Edition},\n"
        "  year      = {2007},\n"
        "  isbn      = {9780131482043},\n"
        "  publisher = {Prentice Hall Press},\n"
        "}\n",
        1 )
  {
    real_type P    = 200;
    real_type Pc   = 33.5;
    real_type T    = 631 * 2;
    real_type Tc   = 126.2;
    real_type Pr   = P / Pc;
    real_type Tr   = T / Tc;
    real_type Asqr = 0.4278 * Pr / pow( Tr, 2.5 );
    real_type B    = 0.0867 * Pr / Tr;
    Q              = B * B + B - Asqr;
    r              = Asqr * B;
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type z = x( 0 );
    f( 0 )      = ( ( z - 1 ) * z - Q ) * z - r;
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    real_type z = x( 0 );
    J.resize( n, n );
    J.insert( 0, 0 ) = ( 3 * z - 2 ) * z - Q;
    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 4 );
    x_vec[0].resize( 1 );
    x_vec[1].resize( 1 );
    x_vec[2].resize( 1 );
    x_vec[3].resize( 1 );
    x_vec[0] << 0.65;
    x_vec[1] << -0.5;
    x_vec[2] << 1;
    x_vec[3] << -0.02;
  }
};
