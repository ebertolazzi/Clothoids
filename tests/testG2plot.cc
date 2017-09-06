#define _USE_MATH_DEFINES
#include "Clothoid.hh"
#include "ClothoidAsyPlot.hh"
#include <cmath>
#include <iostream>
#include <stack>
#include <ctime>

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

int
main() {
	tic();

	Clothoid::G2solve3arc  g2sol;
	
  #if 1
	valueType x0 = -2;
	valueType y0 = 3;
	valueType th0 = M_PI/3;
	valueType k0 = 10.2;
	valueType s0 = M_PI/2/10.2;
	s0 = 0.2;
	valueType x1 = 3;
	valueType y1 = 2;
	valueType th1 = M_PI/10;
	valueType k1 = -0.5;
	valueType s1 = 1;
	#else
	valueType x0  = -1;
	valueType y0  = 0;
	valueType th0 = M_PI / 2;
	valueType k0  = -1;
	valueType s0  = 0.2;
	valueType x1  = 1;
	valueType y1  = 0;
	valueType th1 = -M_PI / 2;
	valueType k1  = -1;
	valueType s1  = 2;
  #endif

	//g2sol.setup(x0, y0, th0, k0, s0, x1, y1, th1, k1, s1);
	g2sol.setup(x0, y0, th0, k0, x1, y1, th1, k1 );
	//g2solve3arc.setup( x0, y0, th0, k0, x1, y1, th1, k1 ) ;
	int iter = g2sol.solve();
	std::cout << "iter = " << iter << '\n' ;

	Clothoid::ClothoidCurve const & S0 = g2sol.getS0();
	Clothoid::ClothoidCurve const & S1 = g2sol.getS1();
	Clothoid::ClothoidCurve const & SM = g2sol.getSM();

	/*
	std::cout << "theta = " << S0.theta(S0.getSmax()) - SM.theta(SM.getSmin()) << '\n' ;
	std::cout << "theta = " << SM.theta(SM.getSmax()) - S1.theta(S1.getSmin()) << '\n' ;

	std::cout << "theta_D = " << S0.theta_D(S0.getSmax()) - SM.theta_D(SM.getSmin()) << '\n' ;
	std::cout << "theta_D = " << SM.theta_D(SM.getSmax()) - S1.theta_D(S1.getSmin()) << '\n' ;

	std::cout << "x = " << S0.X(S0.getSmax()) - SM.X(SM.getSmin()) << '\n' ;
	std::cout << "y = " << SM.Y(SM.getSmax()) - S1.Y(S1.getSmin()) << '\n' ;
	*/

	toc();
	std::cout << "\n\nS0 (NEW)\n" << S0;
	std::cout << "\n\nSM (NEW)\n" << SM;
	std::cout << "\n\nS1 (NEW)\n" << S1;
	//std::cout << "kappa1 = " << S0.getKappa() + S0.getKappa_D()*S0.getL() << '\n';
	//std::cout << "kappa2 = " << S1.getKappa() << '\n';
	
	{
		Clothoid::AsyPlot ap("G2Interpolation.asy", false);
		ap.drawClothoid(S0, "red+bp");
		ap.drawClothoid(SM, "blue+bp");
		ap.drawClothoid(S1, "red+bp");
		valueType factor = 1.1;
		ap.displayAxes("x", "y", x0*factor, x1*factor, y0*factor, y1*factor);
	}
	

	system("pause");

	return 0;
}
