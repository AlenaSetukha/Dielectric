#include <iostream>
#include <cmath>
#include <complex>

#include "common_type.h"
#include "constants.h"
#include "element_geom.h"

//=====================================================Набор правой части==================================================
void get_b(std::complex<double>* b, double* e0, double* k_vec, const TGrid_DC_Full& a)
{
    Constants c;
    int i0;
    std::complex<double> deg, deg1, e_inc[3];

    for (int i = 0; i < a.num_frm; i++)
    {
        i0 = 2 * i;
        deg = scal_prod(k_vec, a.rkt[i]) * c.i_complex;
        deg1 = pow(c.e, deg);

        e_inc[0] = e0[0] * deg1;
        e_inc[1] = e0[1] * deg1;
        e_inc[2] = e0[2] * deg1;
        b[i0] = -scal_prod(e_inc, a.tau[i][0]);
        b[i0 + 1] = -scal_prod(e_inc, a.tau[i][1]);
    }
    return;
}

