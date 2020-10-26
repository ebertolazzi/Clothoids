//#define _USE_MATH_DEFINES
#include "Clothoids.hh"
#include <stack>
#include <ctime>

using G2lib::real_type;
using namespace std;

int
main() {

  static const real_type m_pi = 3.14159265358979323846264338328;

  map<int,int> stats;

  G2lib::G2solve3arc g2solve3arc;
  Utils::TicToc      tictoc;

  int NMAX = 128/8;
  int nkur = 64/8;

  real_type x0 = 0;
  real_type y0 = 0;
  real_type x1 = 1;
  real_type y1 = 0;

  // insert code here...
  real_type thmin = -m_pi*0.999;
  real_type thmax =  m_pi*0.999;

  real_type kur[1000], kmax = 10;
  real_type a = exp( 2*log(kmax)/(nkur-1) );
  nkur = 2*nkur+1;
  cout << "a = " << a << "\n";
  real_type k0 = 1/kmax;
  kur[0] = 0;
  for ( int ii = 1; ii < nkur; ii += 2 ) {
    kur[ii]   = k0;
    kur[ii+1] = -kur[ii];
    k0 *= a;
  }
  cout << "\nkur[0]      = " << kur[0]
       << "\nkur[1]      = " << kur[1]
       << "\nkur[2]      = " << kur[2]
       << "\nkur[nkur-1] = " << kur[nkur-1]
       << "\nkur[nkur-2] = " << kur[nkur-2]
       << "\n";

  // real_type kur[] = {-1e3, -100,-10,-1,-0.1,-0.01,-0.001,-0.0001,0, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1e3 };
  tictoc.tic();
  for ( int ii = 0; ii < nkur; ++ii ) {
    cout << "ii = " << ii << '\n';
    k0 = kur[ii];
    for ( int jj = 0; jj < nkur; ++jj ) {
      real_type k1 = kur[jj];
      for ( int i = 0; i < NMAX; ++i ) {
        real_type th0 = thmin + ((thmax-thmin)*i)/(NMAX-1);
        for ( int j = 0; j < NMAX; ++j ) {
          real_type th1 = thmin + ((thmax-thmin)*j)/(NMAX-1);
          int iter = g2solve3arc.build( x0, y0, th0, k0, x1, y1, th1, k1 );
          if ( iter < 0 ) {
            cout << "iter = " << iter << '\n';
            iter = g2solve3arc.build( x0, y0, th0, k0, x1, y1, th1, k1 );
            cout << "iter dopo = " << iter << '\n';
          }
          stats[iter]++;
        }
      }
    }
  }
  tictoc.toc();
  cout << "stats\n";
  long N = 0;
  for ( map<int,int>::const_iterator is = stats.begin();
        is != stats.end(); ++is ) {
    cout << "iter = " << is->first << " -- " << is->second << '\n';
    N += is->second;
  }
  cout
   << "elapsed = " << tictoc.elapsed_s() << "[s]\n"
   << "ave     = " << tictoc.elapsed_s()/N << "[s]\n"
   << "All done\n";
  return 0;
}

