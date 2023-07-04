#ifndef _INTEGRAL_UNIVERSAL4_ABS_H_
#define _INTEGRAL_INIVERSAL4_ABS_H_


#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>

#include "f_par.h"
#include "integral_par.h"
#include "element_geom.h"


//---------------Description-----------------
// x - точка коллокации
// rut0 - ячейка(4 вершины)
// f_0 - функция ядра
// param - параметры функции(k, радиус сглаживания, машинный ноль, доп параметры )
// int_param - параметры вычисления интеграла(размерность, стартовое разбиение, предельное разбиенеие, точность)
// res - результат вычисления интеграла
//-------------------------------------------

template <typename P>
void integral_universal4_abs(double* x, double** rut0,
        void (*f_0)(double*, double*, f_par, P*), f_par param,
        integral_par int_param, P* res)
{
    double p, q, p1, q1, s, a[3], b[3], a1[3], a2[3], a3[3], a4[3], m1[3], m2[3], rc[3], rn[3];
    double delta = 0.;
    int n = int_param.n_start;
    int p_n;

    P* ff = new P[int_param.idim];
    P* res_prev = new P[int_param.idim];
    for (int g = 0; g < int_param.idim; g++)
    {
        ff[g] = static_cast<P>(0);
        res[g] = static_cast<P>(0);
        res_prev[g] = static_cast<P>(0);
    }


    for (p_n = 0; p_n < int_param.p_max; p_n++)
    {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                p = (double)i / (double)n;
                q = (double)j / (double)n;
                p1 = ((double)i + 1.) / (double)n;
                q1 = ((double)j + 1.) / (double)n;

                for (int k = 0; k < 3; k++) {
                    a[k] = q * rut0[1][k] + (1 - q) * rut0[0][k];
                    b[k] = q * rut0[2][k] + (1 - q) * rut0[3][k];
                    a1[k] = p * b[k] + (1 - p) * a[k];
                    a4[k] = p1 * b[k] + (1 - p1) * a[k];
                }

                for (int k = 0; k < 3; k++) {
                    a[k] = q1 * rut0[1][k] + (1 - q1) * rut0[0][k];
                    b[k] = q1 * rut0[2][k] + (1 - q1) * rut0[3][k];
                    a2[k] = p * b[k] + (1 - p) * a[k];
                    a3[k] = p1 * b[k] + (1 - p1) * a[k];
                }

                for (int k = 0; k < 3; k++) {
                    rc[k] = (a1[k] + a2[k] + a3[k] + a4[k]) / 4.0;
                    m1[k] = ((a2[k] + a3[k]) - (a1[k] + a4[k])) / 2.0;
                    m2[k] = ((a3[k] + a4[k]) - (a1[k] + a2[k])) / 2.0;
                }

                vec_prod(m1, m2, rn);
                s = vec_length(rn);

                f_0(x, rc, param, ff);

                for (int g = 0; g < int_param.idim; g++)
                {
                    res[g] += ff[g] * static_cast<P>(s);
                }
            }
        }

        delta = 0.;
        for (int g = 0; g < int_param.idim; g++)
        {
            delta += abs_tmp(res[g] - res_prev[g]) * abs_tmp(res[g] - res_prev[g]);
        }

        if (delta <  int_param.eps * int_param.eps && p_n != 0)
        {
            break;
        }

        n = n * 2;//шаг уменьшается в 2 раза
        for (int g = 0; g < int_param.idim; g++)
        {
            res_prev[g] = res[g];
            res[g] = static_cast<P>(0);
        }
    }

    if (p_n == int_param.p_max)
    {
        for (int g = 0; g < int_param.idim; g++)
        {
            res[g] = res_prev[g];
        }
    }



    delete[] res_prev;
    delete[] ff;
    return;
}
#endif
