#ifndef _F_PAR_H_
#define _F_PAR_H_

#include <complex>

class f_par
{
public:
    double rs; //радиус сглаживания
    std::complex<double> k;// волновое число
    double n[3];//вектор внешней нормали(если есть/нужен)
    double a[3]; // вектор(a, tau,...)

    f_par(double rs_in, std::complex<double> k_in);

    f_par(double rs_in, std::complex<double> k_in, double* n_in);
};
#endif
