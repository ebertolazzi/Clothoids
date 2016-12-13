#define _USE_MATH_DEFINES
#include "Clothoid.hh"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using Clothoid::valueType ;

int
main( int argc, char const * argv [] ) {

  Clothoid::G2solve3arc g2solve3arc0, g2solve3arc1, g2solve3arc2, g2solve3arc3 ;

  // valori fissi
  valueType x0  = -1 ;
  valueType y0  =  0 ;
  valueType x1  =  1 ;
  valueType y1  =  0 ;

  valueType th_min     = -M_PI ;
  valueType th_max     =  M_PI ;
  valueType k_exp_min  = -4 ;
  valueType k_exp_max  =  2 ;

  // dati per tabella piccola
  int N_TH = 4 ; //35 ;
  int N_K  = 4 ; //20 ;
  
  std::ofstream file0("dataPP.txt") ;
  std::ofstream file1("dataMP.txt") ;
  std::ofstream file2("dataPM.txt") ;
  std::ofstream file3("dataMM.txt") ;
  file0 << "theta0\ttheta1\tkappa0\tkappa1\talpha\tbeta\n" ;
  file1 << "theta0\ttheta1\tkappa0\tkappa1\talpha\tbeta\n" ;
  file2 << "theta0\ttheta1\tkappa0\tkappa1\talpha\tbeta\n" ;
  file3 << "theta0\ttheta1\tkappa0\tkappa1\talpha\tbeta\n" ;
  file0.precision(20) ;
  file1.precision(20) ;
  file2.precision(20) ;
  file3.precision(20) ;

  for ( int i_th0 = 1 ; i_th0 < N_TH ; ++i_th0 ) {
    valueType sth0 = i_th0/valueType(N_TH+1) ;
    valueType th0  = th_min*(1-sth0) + th_max*sth0 ;
    std::cout << "th0 = " << th0 << "\n" ;
    for ( int i_th1 = 1 ; i_th1 < N_TH ; ++i_th1 ) {
      valueType sth1 = i_th1/valueType(N_TH+1) ;
      valueType th1  = th_min*(1-sth1) + th_max*sth1 ;
      std::cout << "done = " << sth0*100 << "% [" << i_th1 << "]\n" ;
      for ( int i_k0 = 0 ; i_k0 < N_K ; ++i_k0 ) {
        valueType sk0 = i_k0/valueType(N_K-1) ;
        valueType k0  = exp( k_exp_min*(1-sk0) + k_exp_max*sk0 ) ;
        for ( int i_k1 = 0 ; i_k1 < N_K ; ++i_k1 ) {
          valueType sk1 = i_k1/valueType(N_K-1) ;
          valueType k1  = exp( k_exp_min*(1-sk1) + k_exp_max*sk1 ) ;
          valueType target0[7], alpha0[7], beta0[7] ;
          valueType target1[7], alpha1[7], beta1[7] ;
          valueType target2[7], alpha2[7], beta2[7] ;
          valueType target3[7], alpha3[7], beta3[7] ;
          bool ok0 = g2solve3arc0.optimize( x0, y0, th0,  k0, x1, y1, th1,  k1, target0, alpha0, beta0 ) ;
          bool ok1 = g2solve3arc1.optimize( x0, y0, th0, -k0, x1, y1, th1,  k1, target1, alpha1, beta1 ) ;
          bool ok2 = g2solve3arc2.optimize( x0, y0, th0,  k0, x1, y1, th1, -k1, target2, alpha2, beta2 ) ;
          bool ok3 = g2solve3arc3.optimize( x0, y0, th0, -k0, x1, y1, th1, -k1, target3, alpha3, beta3 ) ;
          if ( ! ( ok0 && ok1 && ok2 && ok3 ) ) {
            std::cerr
              << "Failed:"
              << "\nth0 = " << th0
              << "\nth1 = " << th1
              << "\nk0  = " << k0
              << "\nk1  = " << k1
              << "\nok0 = " << (ok0?"TRUE":"FALSE")
              << "\nok1 = " << (ok1?"TRUE":"FALSE")
              << "\nok2 = " << (ok2?"TRUE":"FALSE")
              << "\nok3 = " << (ok3?"TRUE":"FALSE")
              << "\n" ;
          } else {
            int kkk = 3 ;
            std::cout
              << th0 << '\t' << th1  << '\t'
              << k0  << '\t' << k1   << '\t'
              << alpha0[kkk] << '\t' << beta0[kkk] << '\n' ;
            file0
              << th0 << '\t' << th1 << '\t'
              << k0  << '\t' << k1  << '\t'
              << alpha0[kkk] << '\t' << beta0[kkk] << '\n' ;
            file1
              << th0 << '\t' << th1 << '\t'
              << -k0 << '\t' << k1  << '\t'
              << alpha1[kkk] << '\t' << beta1[kkk] << '\n' ;
            file2
              << th0 << '\t' << th1 << '\t'
              << k0  << '\t' << -k1 << '\t'
              << alpha2[kkk] << '\t' << beta2[kkk] << '\n' ;
            file3
              << th0 << '\t' << th1 << '\t'
              << -k0 << '\t' << -k1 << '\t'
              << alpha3[kkk] << '\t' << beta3[kkk] << '\n' ;
          }
        }
      }
    }
  }
  file0.close() ;
  file1.close() ;
  file2.close() ;
  file3.close() ;
  return 0 ;
}
