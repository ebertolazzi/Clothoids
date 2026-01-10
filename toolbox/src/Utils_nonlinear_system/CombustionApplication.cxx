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

class CombustionApplication : public NonlinearSystem
{
public:
  // sum log(xi-2)^2+log(xi-10)^2 - prod( xi) ^(1/5)
  CombustionApplication()
    : NonlinearSystem(
        "Combustion Application",
        "@article{Grosan:2012,\n"
        "  title   = {SOLVING POLYNOMIAL SYSTEMS USING A MODIFIED LINE "
        "SEARCH APPROACH},\n"
        "  author  = {Crina Grosan and Ajith Abraham and Vaclav Snasel},\n"
        "  journal = {International Journal of Innovative Computing, "
        "Information and Control},\n"
        "  volume  = {8},\n"
        "  number  = {1},\n"
        "  year    = {2012}\n"
        "}\n\n"
        "@book{Morgan:2009,\n"
        "  author = {Morgan, A.},\n"
        "  title  = {Solving Polynomial Systems Using Continuation for\n"
        "            Engineering and Scientific Problems},\n"
        "  publisher = {Society for Industrial and Applied Mathematics},\n"
        "  year = {2009},\n"
        "  doi = {10.1137/1.9780898719031},\n"
        "}\n\n"
        "@article{Hentenryck:1997,\n"
        "  author  = {Van Hentenryck, P. and McAllester, D. and Kapur, "
        "D.},\n"
        "  title   = {Solving Polynomial Systems Using a Branch and Prune "
        "Approach},\n"
        "  journal = {SIAM Journal on Numerical Analysis},\n"
        "  year    = {1997},\n"
        "  volume  = {34},\n"
        "  number  = {2},\n"
        "  pages   = {797-827},\n"
        "  doi     = {10.1137/S0036142995281504}\n"
        "}\n",
        10 )
  {
  }

  virtual void evaluate( Vector const & x, Vector & f ) const override
  {
    f( 0 ) = x( 1 ) + 2 * x( 5 ) + x( 8 ) + 2 * x( 9 ) - 1e-5;
    f( 1 ) = x( 2 ) + x( 7 ) - 3e-5;
    f( 2 ) = x( 0 ) + x( 2 ) + 2 * x( 4 ) + 2 * x( 7 ) + x( 8 ) + x( 9 ) - 5e-5;
    f( 3 ) = x( 3 ) + 2 * x( 6 ) - 1e-5;
    f( 4 ) = 0.5140437e-7 * x( 4 ) - x( 0 ) * x( 0 );
    f( 5 ) = 0.1006932e-6 * x( 5 ) - 2 * x( 1 ) * x( 1 );
    f( 6 ) = 0.7816278e-15 * x( 6 ) - x( 3 ) * x( 3 );
    f( 7 ) = 0.1496236e-6 * x( 7 ) - x( 0 ) * x( 2 );
    f( 8 ) = 0.6194411e-7 * x( 8 ) - x( 0 ) * x( 1 );
    f( 9 ) = 0.2089296e-14 * x( 9 ) - x( 0 ) * x( 1 );
  }

  virtual void jacobian( Vector const & x, SparseMatrix & J ) const override
  {
    J.resize( n, n );
    J.setZero();

    J.insert( 0, 1 ) = 1;
    J.insert( 0, 5 ) = 2;
    J.insert( 0, 8 ) = 1;
    J.insert( 0, 9 ) = 2;

    J.insert( 1, 2 ) = 1;
    J.insert( 1, 7 ) = 1;

    J.insert( 2, 0 ) = 1;
    J.insert( 2, 2 ) = 1;
    J.insert( 2, 4 ) = 2;
    J.insert( 2, 7 ) = 2;
    J.insert( 2, 8 ) = 1;
    J.insert( 2, 9 ) = 1;

    J.insert( 3, 3 ) = 1;
    J.insert( 3, 6 ) = 2;

    J.insert( 4, 4 ) = 0.5140437e-7;
    J.insert( 4, 0 ) = -2 * x( 0 );

    J.insert( 5, 5 ) = 0.1006932e-6;
    J.insert( 5, 1 ) = -4 * x( 1 );

    J.insert( 6, 6 ) = 0.7816278e-15;
    J.insert( 6, 3 ) = -2 * x( 3 );

    J.insert( 7, 7 ) = 0.1496236e-6;
    J.insert( 7, 0 ) = -x( 2 );
    J.insert( 7, 2 ) = -x( 0 );

    J.insert( 8, 8 ) = 0.6194411e-7;
    J.insert( 8, 0 ) = -x( 1 );
    J.insert( 8, 1 ) = -x( 0 );

    J.insert( 9, 9 ) = 0.2089296e-14;
    J.insert( 9, 0 ) = -x( 1 );
    J.insert( 9, 1 ) = -x( 0 );

    J.makeCompressed();
  }

  virtual void initial_points( vector<Vector> & x_vec ) const override
  {
    x_vec.resize( 1 );
    auto & x0{ x_vec[0] };
    x0.resize( n );
    x0.fill( 1 );
  }
};
