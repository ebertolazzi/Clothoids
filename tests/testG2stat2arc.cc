#define _USE_MATH_DEFINES
#include "Clothoid.hh"
#include "ClothoidAsyPlot.hh"
#include "TicToc.hh"
#include <cmath>
#include <iostream>
#include <stack>
#include <ctime>
#include <map>

std::stack<clock_t> tictoc_stack;

void tic() {
	tictoc_stack.push(clock());
}

void toc() {
	std::cout << "Time elapsed: "
		<< ((double)(clock() - tictoc_stack.top())) / CLOCKS_PER_SEC
		<< std::endl;
	tictoc_stack.pop();
}

using Clothoid::valueType;

static const valueType m_pi = 3.14159265358979323846264338328;

using namespace std ;

map<int,int> stats ;

int
main(int argc, const char * argv[]) {

  Clothoid::G2solve2arc g2solve2arc ;
  TicToc                tictoc ;

  int NMAX = 128/8 ;

  valueType x0 = 0 ;
  valueType y0 = 0 ;
  valueType x1 = 1 ;
  valueType y1 = 0 ;

  // insert code here...
  valueType thmin = -m_pi*0.99 ;
  valueType thmax =  m_pi*0.99 ;

  valueType kur[1000], kmax = 10 ;
  int nkur = 64 ;
  valueType a = exp( 2*log(kmax)/(nkur-1) ) ;
  nkur = 2*nkur+1 ;
  cout << "a = " << a << "\n" ;
  valueType k0 = 1/kmax ;
  kur[0] = 0 ;
  for ( int ii = 1 ; ii < nkur ; ii += 2 ) {
    kur[ii]   = k0 ;
    kur[ii+1] = -kur[ii] ;
    k0 *= a;
  }

  cout << "\nkur[0]      = " << kur[0]
       << "\nkur[1]      = " << kur[1]
       << "\nkur[2]      = " << kur[2]
       << "\nkur[nkur-1] = " << kur[nkur-1]
       << "\nkur[nkur-2] = " << kur[nkur-2]
       << "\n" ;


  // valueType kur[] = {-1e3, -100,-10,-1,-0.1,-0.01,-0.001,-0.0001,0, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1e3 } ;
  tictoc.tic();
  for ( int ii = 0 ; ii < nkur ; ++ii ) {
    valueType k0 = kur[ii];
    for ( int jj = 0 ; jj < nkur ; ++jj ) {
      valueType k1 = kur[jj] ;
      for ( int i = 0 ; i < NMAX ; ++i ) {
        valueType th0 = thmin + ((thmax-thmin)*i)/(NMAX-1);
        for ( int j = 0 ; j < NMAX ; ++j ) {
          valueType th1 = thmin + ((thmax-thmin)*j)/(NMAX-1);
          int iter = g2solve2arc.build( x0, y0, th0, k0, x1, y1, th1, k1 ) ;
          //if ( iter < 0 ) {
          //  cout << "iter = " << iter << '\n' ;
          // iter = g2solve2arc.build( x0, y0, th0, k0, x1, y1, th1, k1 ) ;
          //  cout << "iter dopo = " << iter << '\n' ;
          //}
          stats[iter]++ ;
        }
      }
    }
  }
  tictoc.tic();

  cout << "stats\n" ;
  int acc = 0, acc1 = 0 ;
  for ( map<int,int>::const_iterator is = stats.begin();
        is != stats.end() ; ++is ) {
    cout << "iter = " << is->first << " -- " << is->second << '\n' ;
    if ( is->first>0 ) acc  += is->second ;
    else               acc1 += is->second ;
  }
  cout << "ok = " << acc << " perc = " << double(acc)/(acc+acc1) << '\n' ;
  cout << "elapsed = " << tictoc.totalElapsedSeconds() << '\n' ;
  cout << "All done\n" ;
  return 0;
}
