//#define _USE_MATH_DEFINES
#include "Clothoids.hh"
#include "Clothoids_fmt.hh"

using G2lib::real_type;
using namespace std;
static constexpr real_type m_pi = 3.14159265358979323846264338328;

int
main() {

  map<int,int> stats;

  G2lib::G2solveCLC g2solveCLC;
  Utils::TicToc     tictoc;

  int const NMAX { 32/8 };
  int       nkur { 32/8 };

  real_type x0 = 0;
  real_type y0 = 0;
  real_type x1 = 1;
  real_type y1 = 0;

  // insert code here...
  real_type const thmin = -m_pi*0.99;
  real_type const thmax{ m_pi*0.99 };

  real_type       kur[1000];
  real_type const kmax{10};
  real_type       a{ exp( 2*log(kmax)/(nkur-1) ) };
  nkur = 2*nkur+1;
  fmt::print( "a = {}\n", a );
  real_type k0 = 1/kmax;
  kur[0] = 0;
  for ( int ii = 1; ii < nkur; ii += 2 ) {
    kur[ii]   = k0;
    kur[ii+1] = -kur[ii];
    k0 *= a;
  }

  fmt::print(
    "kur[0]      = {}\n"
    "kur[1]      = {}\n"
    "kur[2]      = {}\n"
    "kur[nkur-1] = {}\n"
    "kur[nkur-2] = {}\n",
    kur[0], kur[1], kur[2], kur[nkur-1], kur[nkur-2]
  );

  // real_type kur[] = {-1e3, -100,-10,-1,-0.1,-0.01,-0.001,-0.0001,0, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1e3 };
  tictoc.tic();
  for ( int ii = 0; ii < nkur; ++ii ) {
    fmt::print( "ii = {}\n", ii );
    k0 = kur[ii];
    for ( int jj{0}; jj < nkur; ++jj ) {
      real_type const k1{ kur[jj] };
      for ( int i{0}; i < NMAX; ++i ) {
        real_type const th0{ thmin + ((thmax-thmin)*i)/(NMAX-1) };
        for ( int j{0}; j < NMAX; ++j ) {
          real_type const th1{ thmin + ((thmax-thmin)*j)/(NMAX-1) };
          int iter{ g2solveCLC.build( x0, y0, th0, k0, x1, y1, th1, k1 ) };
          //if ( iter < 0 ) {
          //  cout << "iter = " << iter << '\n';
          // iter = g2solve2arc.build( x0, y0, th0, k0, x1, y1, th1, k1 );
          //  cout << "iter dopo = " << iter << '\n';
          //}
          stats[iter]++;
        }
      }
    }
  }
  tictoc.tic();

  fmt::print( "stats\n" );
  int acc = 0, acc1 = 0;
  for ( const auto &[fst, snd] : stats ) {
    fmt::print( "iter = {} -- {}\n", fst, snd );
    if ( fst>0 ) acc  += snd;
    else         acc1 += snd;
  }
  fmt::print(
    "ok = {} perc = {}\n"
    "elapsed = {}\n"
    "All done\n",
    acc, static_cast<double>(acc)/(acc+acc1), tictoc.elapsed_s()
  );
  return 0;
}
