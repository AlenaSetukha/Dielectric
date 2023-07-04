#ifndef _INTEGRAL_PAR_H_
#define _INTEGRAL_PAR_H_

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>

class integral_par
{
public:
    int idim;// размерность функции
    int n_start;// стартовое разбиение ячейки
    int p_max;// предельное число разбиений(итераций)- 2^p_max
    double eps;// точность вычисления интеграла(имеет смысл при p_max != 1)

    integral_par(int idim_in, int n_start_in, int p_max_in, double eps_in);

    integral_par(const integral_par& obj);
};
#endif

