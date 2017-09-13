#define _USE_MATH_DEFINES
#include "Clothoid.hh"
#include "ClothoidAsyPlot.hh"
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

  Clothoid::G2solve3arc g2solve3arc ;

  int NMAX = 32 ;

  valueType x0 = 0 ;
  valueType y0 = 0 ;
  valueType x1 = 1 ;
  valueType y1 = 0 ;

  // insert code here...
  valueType thmin = -m_pi*0.99 ;
  valueType thmax =  m_pi*0.99 ;

  int nkur = 17 ;
  valueType kur[] = {-1e3, -100,-10,-1,-0.1,-0.01,-0.001,-0.0001,0, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1e3 } ;
  for ( int ii = 0 ; ii < nkur ; ++ii ) {
    valueType k0 = kur[ii];
    for ( int jj = 0 ; jj < nkur ; ++jj ) {
      valueType k1 = kur[jj] ;
      for ( int i = 0 ; i < NMAX ; ++i ) {
        valueType th0 = thmin + ((thmax-thmin)*i)/(NMAX-1);
        for ( int j = 0 ; j < NMAX ; ++j ) {
          valueType th1 = thmin + ((thmax-thmin)*j)/(NMAX-1);
          int iter = g2solve3arc.build( x0, y0, th0, k0, x1, y1, th1, k1 ) ;
          if ( iter < 0 ) {
            cout << "iter = " << iter << '\n' ;
            iter = g2solve3arc.build( x0, y0, th0, k0, x1, y1, th1, k1 ) ;
            cout << "iter dopo = " << iter << '\n' ;
          }
          stats[iter]++ ;
        }
      }
    }
  }
  cout << "stats\n" ;
  for ( map<int,int>::const_iterator is = stats.begin();
        is != stats.end() ; ++is )
    cout << "iter = " << is->first << " -- " << is->second << '\n' ;
  return 0;
}
