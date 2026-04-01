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

class Beale : public NonlinearSystem
{
public:
  Beale()
    : NonlinearSystem(
        "Beale",
        "@book{beale1958,\n"
        "  title    = {On an Iterative Method for Finding a Local Minimum\n"
        "              of a Function of More Than One Variable},\n"
        "  author    = {Beale, E.M.L.},\n"
        "  series    = {Technical report\n"
        "               (Princeton University. Statistical Techniques "
        "Research Group)},\n"
        "  year      = {1958},\n"
        "  publisher = {Statistical Techniques Research Group,\n"
        "               Section of Mathematical Statistics,\n"
        "               Department of Mathematics, Princeton University}\n"
        "}\n\n"
        "@book{brent2013,\n"
        "  author    = {Brent, R.P.},\n"
        "  title     = {Algorithms for Minimization Without Derivatives},\n"
        "  isbn      = {9780486143682},\n"
        "  series    = {Dover Books on Mathematics},\n"
        "  year      = {2013},\n"
        "  publisher = {Dover Publications}\n"
        "}\n",
        2 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    real_type x1 = x( 0 );
    real_type x2 = x( 1 );

    real_type f1 = 1.5 - x1 * ( 1.0 - x2 );
    real_type f2 = 2.25 - x1 * ( 1.0 - x2 * x2 );
    real_type f3 = 2.625 - x1 * ( 1.0 - x2 * x2 * x2 );

    real_type df1dx1 = x2 - 1;
    real_type df1dx2 = x1;
    real_type df2dx1 = x2 * x2 - 1;
    real_type df2dx2 = 2.0 * x1 * x2;
    real_type df3dx1 = x2 * x2 * x2 - 1.0;
    real_type df3dx2 = 3.0 * x1 * x2 * x2;

    f( 0 ) = 2.0 * ( f1 * df1dx1 + f2 * df2dx1 + f3 * df3dx1 );
    f( 1 ) = 2.0 * ( f1 * df1dx2 + f2 * df2dx2 + f3 * df3dx2 );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    real_type x1 = x( 0 );
    real_type x2 = x( 1 );

    real_type f1 = 1.5 - x1 * ( 1.0 - x2 );
    real_type f2 = 2.25 - x1 * ( 1.0 - x2 * x2 );
    real_type f3 = 2.625 - x1 * ( 1.0 - x2 * x2 * x2 );

    real_type df1dx1 = x2 - 1;
    real_type df1dx2 = x1;
    real_type df2dx1 = x2 * x2 - 1;
    real_type df2dx2 = 2.0 * x1 * x2;
    real_type df3dx1 = x2 * x2 * x2 - 1.0;
    real_type df3dx2 = 3.0 * x1 * x2 * x2;

    real_type d2f1dx12 = 1.0;
    real_type d2f1dx21 = 1.0;

    real_type d2f2dx12 = 2.0 * x2;
    real_type d2f2dx21 = 2.0 * x2;
    real_type d2f2dx22 = 2.0 * x1;

    real_type d2f3dx12 = 3.0 * x2 * x2;
    real_type d2f3dx21 = 3.0 * x2 * x2;
    real_type d2f3dx22 = 6.0 * x1 * x2;

    J.insert( 0, 0 ) = 2.0 * ( df1dx1 * df1dx1 + df2dx1 * df2dx1 + df3dx1 * df3dx1 );

    J.insert( 0, 1 ) = 2.0 * ( df1dx2 * df1dx1 + f1 * d2f1dx12 + df2dx2 * df2dx1 + f2 * d2f2dx12 + df3dx2 * df3dx1 +
                               f3 * d2f3dx12 );

    J.insert( 1, 0 ) = 2.0 * ( df1dx1 * df1dx2 + f1 * d2f1dx21 + df2dx1 * df2dx2 + f2 * d2f2dx21 + df3dx1 * df3dx2 +
                               f3 * d2f3dx21 );

    J.insert( 1, 1 ) = 2.0 * ( df1dx2 * df1dx2 + df2dx2 * df2dx2 + f2 * d2f2dx22 + df3dx2 * df3dx2 + f3 * d2f3dx22 );
    J.makeCompressed();
  }

  virtual void exact_solution( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 3, 0.5;
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0 << 1, 1;
  }
};
