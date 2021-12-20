//#define _USE_MATH_DEFINES
#include "Clothoids.hh"

using G2lib::real_type;
using G2lib::int_type;
using namespace std;

int
main() {

  // check constructors
  real_type x0     = 0;
  real_type y0     = 2;
  real_type x1     = 4;
  real_type y1     = 3.5;
  real_type x2     = 6;
  real_type y2     = 1;


  G2lib::Triangle2D T1;
  T1.build( x0, y0, x1, y1, x2, y2, 0, 0, 0 );

  int_type icode = T1.isInside( T1.baricenterX(), T1.baricenterY() );

  cout << "icode = " << icode << "\n";

  icode = T1.isInside( T1.P1() );

  cout << "icode = " << icode << "\n";

  icode = T1.isInside( 10, 20 );

  cout << "icode = " << icode << "\n";

  {
    G2lib::ClothoidCurve C;
    real_type xx0    = 0;
    real_type yy0    = 2;
    real_type theta0 = 0;
    real_type kappa0 = 10;
    real_type dk     = -1;// 0;
    real_type L      = 1;
    C.build( xx0, yy0, theta0, kappa0, dk, L );

    vector<G2lib::Triangle2D> tvec;
    C.bbTriangles_ISO( 0, tvec );

    for ( size_t i = 0; i < tvec.size(); ++i )
      cout << i << " " << tvec[i] << "\n\n";

    cout << tvec.size() << '\n';
  }

  cout << "\n\nALL DONE FOLKS!!!\n";

  return 0;
}
