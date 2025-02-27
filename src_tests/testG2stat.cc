//#define _USE_MATH_DEFINES
#include "Clothoids.hh"
#include "Clothoids_fmt.hh"

#include <stack>
#include <ctime>

using G2lib::real_type;
using namespace std;

int
main() {

  static constexpr real_type m_pi = 3.14159265358979323846264338328;

  map<int,int> stats;

  G2lib::G2solve3arc g2solve3arc;
  Utils::TicToc      tictoc;

  constexpr int NMAX { 128/8 };
  int       nkur { 64/8  };

  constexpr real_type x0 { 0 };
  constexpr real_type y0 { 0 };
  constexpr real_type x1 { 1 };
  constexpr real_type y1 { 0 };

  // insert code here...
  real_type const thmin { -m_pi*0.999 };
  real_type const thmax {  m_pi*0.999 };

  real_type           kur[1000];
  constexpr real_type kmax { 10 };
  real_type           a    { exp( 2*log(kmax)/(nkur-1) ) };
  nkur = 2*nkur+1;
  fmt::print( "a = {}\n", a );
  real_type k0{ 1/kmax };
  kur[0] = 0;
  for ( int ii{1}; ii < nkur; ii += 2 ) {
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
  for ( int ii{0}; ii < nkur; ++ii ) {
    fmt::print( "ii = {}\n", ii );
    k0 = kur[ii];
    for ( int jj{0}; jj < nkur; ++jj ) {
      real_type const k1{ kur[jj] };
      for ( int i{0}; i < NMAX; ++i ) {
        real_type const th0{ thmin + ((thmax-thmin)*i)/(NMAX-1) };
        for ( int j{0}; j < NMAX; ++j ) {
          real_type const th1{ thmin + ((thmax-thmin)*j)/(NMAX-1) };
          int iter{ g2solve3arc.build( x0, y0, th0, k0, x1, y1, th1, k1 ) };
          if ( iter < 0 ) {
            fmt::print( "iter = {}\n", iter );
            iter = g2solve3arc.build( x0, y0, th0, k0, x1, y1, th1, k1 );
            fmt::print( "iter dopo = {}\n", iter );
          }
          stats[iter]++;
        }
      }
    }
  }
  tictoc.toc();
  fmt::print( "stats\n" );
  long N{0};
  for ( const auto &[fst, snd] : stats ) {
    fmt::print( "iter = {} -- {}\n", fst, snd );
    N += snd;
  }
  fmt::print(
    "elapsed = {} [s]\n"
    "ave     = {} [s]\n"
    "All done\n",
    tictoc.elapsed_s(), tictoc.elapsed_s()/N
  );
  return 0;
}

