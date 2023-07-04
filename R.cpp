#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>

#include "integral_par.h"
#include "f_par.h"
#include "integral_universal4_abs.h"
#include "kernel_lib.h"
#include "element_geom.h"


void R_rot(std::complex<double>* j, double* x, double** rut0,
        integral_par integral_par_f_grad_simple, f_par param,
        std::complex<double>* res)
{
    std::complex<double> cur_res3[3];
    integral_universal4_abs(x, rut0, f_grad_simple_pot_G, param, integral_par_f_grad_simple, cur_res3);
    vec_prod(cur_res3, j, res);
    return;
}
