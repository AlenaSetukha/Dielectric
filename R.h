#ifndef _R_H_
#define _R_H_

#include <complex>
#include "integral_par.h"
#include "f_par.h"

void R_rot(std::complex<double>* j, double* x, double** rut0,
        integral_par integral_par_f_grad_simple, f_par param,
        std::complex<double>* res);
#endif
