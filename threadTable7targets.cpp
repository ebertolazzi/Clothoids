#define _USE_MATH_DEFINES
#include "Clothoid.hh"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <thread>


using Clothoid::valueType;
using Clothoid::indexType;

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
	indexType const N_files = 8;
	indexType sign0 = 0;
	indexType sign1 = 0;
	std::ofstream files[N_files];

	if (thr_numb == 0) {
		files[0].open("dataPP_1_length.txt");
		files[1].open("dataPP_2_curvature.txt");
		files[2].open("dataPP_3_jerk.txt");
		files[3].open("dataPP_4_snap.txt");
		files[4].open("dataPP_5_TV-angle.txt");
		files[5].open("dataPP_6_TV-curvature.txt");
		files[6].open("dataPP_7_target.txt");
		files[7].open("dataPP_failed.txt");
		sign0 = 1;
		sign1 = 1;
	}
	if (thr_numb == 1) {
		files[0].open("dataMP_1_length.txt");
		files[1].open("dataMP_2_curvature.txt");
		files[2].open("dataMP_3_jerk.txt");
		files[3].open("dataMP_4_snap.txt");
		files[4].open("dataMP_5_TV-angle.txt");
		files[5].open("dataMP_6_TV-curvature.txt");
		files[6].open("dataMP_7_target.txt");
		files[7].open("dataMP_failed.txt");
		sign0 = -1;
		sign1 = 1;
	}
	if (thr_numb == 2) {
		files[0].open("dataPM_1_length.txt");
		files[1].open("dataPM_2_curvature.txt");
		files[2].open("dataPM_3_jerk.txt");
		files[3].open("dataPM_4_snap.txt");
		files[4].open("dataPM_5_TV-angle.txt");
		files[5].open("dataPM_6_TV-curvature.txt");
		files[6].open("dataPM_7_target.txt");
		files[7].open("dataPM_failed.txt");
		sign0 = 1;
		sign1 = -1;
	}
	if (thr_numb == 3) {
		files[0].open("dataMM_1_length.txt");
		files[1].open("dataMM_2_curvature.txt");
		files[2].open("dataMM_3_jerk.txt");
		files[3].open("dataMM_4_snap.txt");
		files[4].open("dataMM_5_TV-angle.txt");
		files[5].open("dataMM_6_TV-curvature.txt");
		files[6].open("dataMM_7_target.txt");
		files[7].open("dataMM_failed.txt");
		sign0 = -1;
		sign1 = -1;
	}

	for (int i = 0; i < N_files; ++i) {
		files[i] << "theta0\ttheta1\tkappa0\tkappa1\talpha\tbeta\n";
		files[i].precision(20);
	}
  
  valueType sth0, th0, sth1, th1, sk0, k0, sk1, k1 ;
	valueType target[7], valpha[7], vbeta[7];
	for (int i_th0 = 1; i_th0 < N_TH; ++i_th0) {
		sth0 = i_th0 / valueType(N_TH + 1);
		th0 = th_min*(1 - sth0) + th_max*sth0;
		std::cout << "th0 = " << th0 << "\n";
		for (int i_th1 = 1; i_th1 < N_TH; ++i_th1) {
      sth1 = i_th1 / valueType(N_TH + 1);
			th1 = th_min*(1 - sth1) + th_max*sth1;
			std::cout << "done = " << sth0 * 100 << "% [" << i_th1 << "]\n";
			for (int i_k0 = 0; i_k0 < N_K; ++i_k0) {
        sk0 = i_k0 / valueType(N_K - 1);
				k0 = exp(k_exp_min*(1 - sk0) + k_exp_max*sk0);
				for (int i_k1 = 0; i_k1 < N_K; ++i_k1) {
					sk1 = i_k1 / valueType(N_K - 1);
					k1 = exp(k_exp_min*(1 - sk1) + k_exp_max*sk1);
					std::cout << "Thread N. " << thr_numb <<
						th0 << " " << th1 << " " << sign0*k0 << " " << sign1*k1 << "\n";

					bool ok = g2solve3arc.optimize(x0, y0, th0, sign0*k0, x1, y1, th1, sign1*k1, target, valpha, vbeta);

					if (!ok) {
            std::cout << "errore\n" ;
            std::cout << x0 << '\n' ;
            std::cout << y0 << '\n' ;
            std::cout << th0 << '\n' ;
            std::cout << th1 << '\n' ;
						files[7]
							<< "  th0 = " << th0
							<< "\nth1 = " << th1
							<< "\nk0  = " << k0
							<< "\nk1  = " << k1
							//<< "\nok  = " << (ok ? "TRUE" : "FALSE")
							<< "\n";
					}
					else {
						for (int i = 0; i < N_files-1; ++i) {
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
	int N_TH = 36;
	int N_K  = 2;
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
