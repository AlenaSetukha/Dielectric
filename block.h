#ifndef _BLOCK_H_
#define _BLOCK_H_

#include <complex>
#include "common_type.h"
#include "integral_par.h"
#include "f_par.h"

//Функции получения блоков 2x2 от операторов K = rotrot, R = rot

void get_block2_K(const TGrid_DC_Full& a, int i, int j,
        integral_par integral_par_f_simple, integral_par integral_par_f_grad_simple,
        f_par param, f_par param_seg,
        std::complex<double> (*block2)[2]);

void get_block2_R(const TGrid_DC_Full& a, int i, int j,
        integral_par int_par_f_K, f_par param_f_grad_R,
        std::complex<double> (*block2)[2], int sign);
#endif
