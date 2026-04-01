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

class LinearFunctionFullRank : public NonlinearSystem
{
public:
  LinearFunctionFullRank()
    : NonlinearSystem(
        "Linear function - full rank",
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
        10 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type sumx = x.array().sum();
    f              = x.array() - ( 2 * sumx / n ) + 1;
  }

  virtual void jacobian( Vector const &, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();
    real_type bf = real_type( 2 ) / n;
    for ( integer i{ 0 }; i < n; ++i )
    {
      for ( integer j{ 0 }; j < n; ++j )
      {
        real_type tmp = -bf;
        if ( i == j ) tmp += 1;
        J.insert( i, j ) = tmp;
      }
    }
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 1 );
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( -1 );
  }
};
