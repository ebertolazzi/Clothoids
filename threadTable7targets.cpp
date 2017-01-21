#define _USE_MATH_DEFINES
#include "Clothoid.hh"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <thread>


using Clothoid::valueType;
using Clothoid::indexType;

static int const N_target = 8 ;

void
threadLaunch( int       thr_numb,
              int       N_TH,
              int       N_K,
              valueType th_min,
              valueType th_max,
              valueType k_exp_min,
              valueType k_exp_max) {

	Clothoid::G2solve3arc g2solve3arc;

	// valori fissi
	valueType x0 = -1;
	valueType y0 = 0;
	valueType x1 = 1;
	valueType y1 = 0;
	indexType const N_files = N_target + 1 ;
	indexType sign0 = 0;
	indexType sign1 = 0;
	std::ofstream files[N_files];
  
  char const *suffix[] = {
		"_1_length.txt",
		"_2_curv.txt",
		"_3_TV-angle.txt",
		"_4_length1.txt",
		"_5_curv1.txt",
		"_6_TV-angle1.txt",
		"_7_curv_TV-angle.txt",
		"_8_curv_jerk.txt",
		"_failed.txt"
  } ;

	if (thr_numb == 0) {
    for ( int i = 0 ; i < N_files ; ++i )
  		files[i].open((std::string("dataPP")+suffix[i]).c_str());
		sign0 = 1;
		sign1 = 1;
	}
	if (thr_numb == 1) {
    for ( int i = 0 ; i < N_files ; ++i )
  		files[i].open((std::string("dataMP")+suffix[i]).c_str());
		sign0 = -1;
		sign1 = 1;
	}
	if (thr_numb == 2) {
    for ( int i = 0 ; i < N_files ; ++i )
  		files[i].open((std::string("dataPM")+suffix[i]).c_str());
		sign0 = 1;
		sign1 = -1;
	}
	if (thr_numb == 3) {
    for ( int i = 0 ; i < N_files ; ++i )
  		files[i].open((std::string("dataMM")+suffix[i]).c_str());
		sign0 = -1;
		sign1 = -1;
	}

	for (int i = 0; i < N_files; ++i) {
		files[i] << "theta0\ttheta1\tkappa0\tkappa1\talpha\tbeta\n";
		files[i].precision(20);
	}
	for (int i_th0 = 1; i_th0 <= N_TH; ++i_th0) {
		valueType sth0 = i_th0 / valueType(N_TH + 1);
		valueType th0 = th_min*(1 - sth0) + th_max*sth0;
		std::cout << "th0 = " << th0 << "\n";
		for (int i_th1 = 1; i_th1 <= N_TH; ++i_th1) {
			valueType sth1 = i_th1 / valueType(N_TH + 1);
			valueType th1 = th_min*(1 - sth1) + th_max*sth1;
			std::cout << "done = " << sth0 * 100 << "% [" << i_th1 << "]\n";
			for (int i_k0 = 0; i_k0 < N_K; ++i_k0) {
				valueType sk0 = i_k0 / valueType(N_K - 1);
				valueType k0 = exp(k_exp_min*(1 - sk0) + k_exp_max*sk0);
				for (int i_k1 = 0; i_k1 < N_K; ++i_k1) {
					valueType sk1 = i_k1 / valueType(N_K - 1);
					valueType k1 = exp(k_exp_min*(1 - sk1) + k_exp_max*sk1);
					valueType target[N_target];
					valueType valpha[N_target];
					valueType vbeta[N_target];
					//std::cout << "Thread N. " << thr_numb <<
					//	th0 << " " << th1 << " " << sign0*k0 << " " << sign1*k1 << "\n";

					bool ok = g2solve3arc.optimize(x0, y0, th0, sign0*k0, x1, y1, th1, sign1*k1, target, valpha, vbeta);

					if (!ok) {
						files[N_target] << th0      << '\t' << th1 << '\t'
							              << sign0*k0 << '\t' << sign1*k1 << '\n';
					}
					else {
						for (int i = 0; i < N_target ; ++i) {
							files[i] << th0 << '\t' << th1 << '\t'
							         << sign0*k0 << '\t' << sign1*k1 << '\t'
                       << valpha[i] << '\t' << vbeta[i] << '\n';
            }
					}
				}
			}
		}
	}

	for (int i = 0; i < N_files; ++i) {
		files[i].close();
	}
}

int
main(int argc, char const * argv[]) {

	valueType th_min = -M_PI;
	valueType th_max = M_PI;
	valueType k_exp_min = -4;
	valueType k_exp_max = 2;

	// dati per tabella piccola 36 e 20
	int N_TH = 18 ; // 72;
	int N_K = 10 ; // 40;
	std::cout << "N_TH = " << N_TH << " N_K = " << N_K << "\n";
	
	
	std::thread threadObj0(threadLaunch, 0, N_TH, N_K, th_min, th_max, k_exp_min, k_exp_max);
	std::thread threadObj1(threadLaunch, 1, N_TH, N_K, th_min, th_max, k_exp_min, k_exp_max);
  std::thread threadObj2(threadLaunch, 2, N_TH, N_K, th_min, th_max, k_exp_min, k_exp_max);
	std::thread threadObj3(threadLaunch, 3, N_TH, N_K, th_min, th_max, k_exp_min, k_exp_max);

	threadObj0.join();
	threadObj1.join();
	threadObj2.join();
	threadObj3.join();
	system("pause");
	return 0;
}
