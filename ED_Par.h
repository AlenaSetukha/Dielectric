#ifndef _ED_PAR_H_
#define _ED_PAR_H_

#include <complex>
#include <string>

class ED_Par
{
public:
    int k_medium;//число областей
    std::complex<double>* eps_d;//если они действ, то просто задать (10., 0.) в файле
    std::complex<double>* m_d;
    std::complex<double>* k;
    std::complex<double>* lambda;

    double k_vec[3];
    double alpha, omega, alpha_start, alpha_end;


    ED_Par(std::string filename);

};
#endif
