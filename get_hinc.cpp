#include <complex>
#include "ED_Par.h"
#include "constants.h"
#include "element_geom.h"

void get_hinc(double* x, double* e0, ED_Par ed_param, std::complex<double>* h_inc)
{
    Constants c;
    std::complex<double>  deg, deg1;
    double vvv[3];
    deg = std::complex<double>(1., 0.) / ed_param.omega / ed_param.m_d[0] / c.m0;
    vec_prod(ed_param.k_vec, e0, vvv);
    deg1 = deg * pow(c.e, c.i_complex * scal_prod(ed_param.k_vec, x));
    h_inc[0] = vvv[0] * deg1;
    h_inc[1] = vvv[1] * deg1;
    h_inc[2] = vvv[2] * deg1;
}
