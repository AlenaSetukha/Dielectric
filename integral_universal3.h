#ifndef _UNTEGRAL_UNIVERSAL3_H_
#define _UNTEGRAL_UNIVERSAL3_H_

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>

#include "f_par.h"
#include "integral_par.h"
#include "element_geom.h"


template <typename P>
void integral_universal3(double* x, double** rut0,
        void* (*f_0)(double*, double*, f_par, P*),
        f_par param, integral_par int_par, P* res)
{
    double p, q, p1, q1, s, a[3], b[3], a1[3], a2[3], a3[3], a4[3], m1[3], m2[3], rc[3], rn[3];
    double delta = 0.;
    int n =  int_param.n_start;


    P* ff = new P[int_param.idim];
    P* res_prev = new P[int_param.idim];
    for (int g = 0; g < int_param.idim; g++)
    {
        ff[g] = static_cast<P>(0);
        res[g] = static_cast<P>(0);
        res_prev[g] = static_cast<P>(0);
    }

    s = tr_square(rut0[0], rut0[1], rut0[2]);
    for (int k = 0; k < 3; k++) {
        m1[k] = rut0[1][k] - rut0[0][k];
        m2[k] = rut0[2][k] - rut0[0][k];
    }


    for (int p = 0; p < int_param.p_max; p++)
    {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                for (int k = 0; k < 3; k++) {
                    a1[k] = rut0[0][k] + (i - 1) * m1[k] + (j - 1) * m2[k];
                    a2[k] = a1[k] + m1[k];
                    a3[k] = a1[k] + m2[k];
                    rc[k] = (a1[k] + a2[k] + a3[k]) / 3.0;
                }

                f_0(x, rc, param, ff);
                for (int g = 0; g < idim; g++)
                {
                    res[g] += ff[g];
                }

                if (j < i - 1) {
                    for (int k = 0; k < 3; k++) {
                        a2[k] = a3[k];
                        a3[k] = a2[k] - m1[k];
                        rc[k] = (a1[k] + a2[k] + a3[k]) / 3.0;
                    }
                    f_0(x, rc, param, ff);
                    for (int g = 0; g < idim; g++)
                    {
                        res[g] += ff[g];
                    }
                }
            }
        }

	    delta = 0.;
        for (int g = 0; g < int_param.idim; g++)
        {
            res[g] = res[g] * s / n / n;
            delta += abs_tmp(res[g] - res_prev[g]) * abs_tmp(res[g] - res_prev[g]);
        }

        if (delta <  int_param.eps * int_param.eps && p_n != 0)//вместо кроня eps^2
        {
            break;
        }
        n = n * 2;
        
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

    delete[] ff;
    delete[] res_prev;
    return;
}
#endif

