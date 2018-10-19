//#define _USE_MATH_DEFINES
#include "Clothoid.hh"
#include "ClothoidAsyPlot.hh"
#include <cmath>
#include <iostream>
#include <stack>
#include <ctime>

#include <stdlib.h>

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

using G2lib::real_type;

int
main() {
	tic();

	G2lib::G2solve2arc g2sol;
	
  #if 1
	real_type x0 = -2;
	real_type y0 = 3;
	real_type th0 = M_PI/3;
	real_type k0 = 10.2;
	real_type s0 = M_PI/2/10.2;
	s0 = 0.2;
	real_type x1 = 3;
	real_type y1 = 2;
	real_type th1 = M_PI/10;
	real_type k1 = -0.5;
	//real_type s1 = 1;
	#else
	real_type x0  = -1;
	real_type y0  = 0;
	real_type th0 = M_PI / 2;
	real_type k0  = -1;
	real_type s0  = 0.2;
	real_type x1  = 1;
	real_type y1  = 0;
	real_type th1 = -M_PI / 2;
	real_type k1  = -1;
	real_type s1  = 2;
  #endif

	//int iter = g2sol.build(x0, y0, th0, k0, s0, x1, y1, th1, k1, s1);
	int iter = g2sol.build(x0, y0, th0, k0, x1, y1, th1, k1 );
	//int iter = g2solve3arc.build( x0, y0, th0, k0, x1, y1, th1, k1 );
	std::cout << "iter = " << iter << '\n';

	G2lib::ClothoidCurve const & S0 = g2sol.getS0();
	G2lib::ClothoidCurve const & S1 = g2sol.getS1();

	toc();
	std::cout << "\n\nS0 (NEW)\n" << S0;

	std::cout << "\n\nS1 (NEW)\n" << S1;

  std::cout << "\n\n\n";

  std::cout << "x1     = " << S0.xEnd()   << "\n";
  std::cout << "x0     = " << S1.xBegin() << "\n\n";

  std::cout << "y1     = " << S0.yEnd()   << "\n";
  std::cout << "y0     = " << S1.yBegin() << "\n\n";

  std::cout << "theta1 = " << S0.thetaEnd()   << "\n";
  std::cout << "theta0 = " << S1.thetaBegin() << "\n\n";

	return 0;
}
