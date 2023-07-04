#ifndef _K_H_
#define _K_H_

#include <complex>
#include "integral_par.h"
#include "f_par.h"

void K_rot_rot(std::complex<double>* j, double* x, double** rut0, double* norm,
        integral_par integral_par_f_simple, integral_par integral_par_f_grad_simple,
        f_par param, f_par param_seg,
        std::complex<double>* res);
#endif
