#ifndef _K0_H_
#define _K0_H_

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>

#include "constants.h"
#include "f_par.h"
#include "integral_par.h"
#include "element_geom.h"
#include "integral_universal_seg_abs.h"




//---------------Description-----------------
// e - вектор, для которого считаем(j)
// x - точка коллокации
// norm - нормаль к ячейке
// rut0 - ячейка(4 вершины)
// f_0 - функция ядра
// param - параметры функции(k, радиус сглаживания, машинный ноль, доп параметры )
// int_param - параметры вычисления интеграла(размерность, стартовое разбиение, предельное разбиенеие, точность)
// res - результат вычисления суммы четырех инетгралов
//-------------------------------------------

template <typename P>
void k0(P* e, double* x, double* norm, double** rut0,
        void (*f_0)(double*, double*, f_par, P*), f_par param,
        integral_par int_par, P* res)
{
    Constants c;
    P* ks = new P[int_par.idim];
    double len_nu;
    double nu[3], diff[3];

    for (int g = 0; g < int_par.idim; g++)
    {
        res[g] = static_cast<P>(0);
        ks[g] = static_cast<P>(0);
    }


    for (int s = 0; s < 4; s++)
    {
//std::cout << x[0] << " " << x[1] << " " << x[2] << std::endl;
        for (int g = 0; g < int_par.idim; g++)
        {
            ks[g] = static_cast<P>(0);
            diff[g] = 0.;
        }

        if (s != 3)
        {

            integral_universal_seg_abs(rut0[s], rut0[s + 1], x, f_0, param, int_par, ks);

            for (int g = 0; g < int_par.idim; g++)
            {
                diff[g] = rut0[s + 1][g] - rut0[s][g];//b-a
            }
        } else {
            integral_universal_seg_abs(rut0[3], rut0[0], x, f_0, param, int_par, ks);
            for (int g = 0; g < int_par.idim; g++)
            {
                diff[g] = rut0[0][g] - rut0[3][g];//b-a
            }
        }

        vec_prod(diff, norm, nu);
        len_nu = vec_length(nu);


        if (len_nu > c.machine_zero)
        {
            for (int g = 0; g < int_par.idim; g++)
            {
                nu[g] /= len_nu;
            }
            for (int g = 0; g < int_par.idim; g++)
            {
                ks[g] *= -scal_prod(e, nu);//ks[3] ab
                res[g] += ks[g];
            }
        }
    }

    delete[] ks;
    return;
}
#endif
