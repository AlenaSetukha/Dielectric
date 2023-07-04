#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>

#include "constants.h"
#include "element_geom.h"
#include "common_type.h"
#include "ED_Par.h"

//------------------------------------Функцияя эпр для идеального проводника-------------------------------------
double f_epr_ideal_conductor(double* tau, const TGrid_DC_Full& a,
        std::complex<double> k, std::complex<double>** j, double* e0)
{
    std::complex<double> res[3];
    std::complex<double> s = 0;
    for (int i = 0; i < 3; i++)
    {
        res[i] = std::complex<double>(0., 0.);
    }

    std::complex<double> e_iktx = std::complex<double>(0., 0.);
    Constants c;
    for (int i = 0; i < a.num_frm; i++)
    {
        e_iktx = exp(-c.i_complex * k * scal_prod(tau, a.rkt[i]));
        e_iktx *= k * k * a.s[i];
        s = scal_prod(j[i], tau);
        for (int g = 0; g < 3; g++)
        {
            res[g] += (j[i][g] - tau[g] * s) * e_iktx * c.pi_reverse;
        }
    }
    return vec_length(res) * vec_length(res) * 4. * c.pi / vec_length(e0) / vec_length(e0);
}



//------------------------------------Функцияя эпр для диэлектрика-------------------------------------
double f_epr_dielectric(double* tau, const TGrid_DC_Full& a,
        std::complex<double>** j_E, std::complex<double>** j_M,
        double* e0, ED_Par ed_param)
{
    std::complex<double> res[3];
    std::complex<double> s = 0;
    std::complex<double> v[3];
    for (int i = 0; i < 3; i++)
    {
        res[i] = std::complex<double>(0., 0.);
    }

    std::complex<double> e_iktx = std::complex<double>(0., 0.);
    Constants c;
    for (int i = 0; i < a.num_frm; i++)
    {
        e_iktx = exp(-c.i_complex * ed_param.k[0] * scal_prod(tau, a.rkt[i]));
        e_iktx *= a.s[i];
        s = scal_prod(j_E[i], tau);
        vec_prod(tau, j_M[i], v);

        for (int g = 0; g < 3; g++)
        {
            res[g] += ((ed_param.k[0] * ed_param.k[0] * c.i_complex  / ed_param.omega / ed_param.eps_d[0] / c.eps0) *
             (j_E[i][g] - tau[g] * s) - (c.i_complex  * ed_param.k[0]) * v[g] ) * e_iktx;
        }
    }
    return vec_length(res) * vec_length(res) * c.pi_reverse / vec_length(e0) / vec_length(e0);
}