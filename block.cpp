#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>

#include "common_type.h"
#include "integral_par.h"
#include "f_par.h"
#include "element_geom.h"
#include "K.h"
#include "R.h"


void get_block2_K(const TGrid_DC_Full& a, int i, int j,
        integral_par integral_par_f_simple, integral_par integral_par_f_grad_simple,
        f_par param, f_par param_seg,
        std::complex<double> (*block2)[2])
{
    std::complex<double> cur_res3[3];

    std::complex<double> j_vec_complex[3];
    j_vec_complex[0] = std::complex<double>(a.tau[j][0][0], 0.);
    j_vec_complex[1] = std::complex<double>(a.tau[j][0][1], 0.);
    j_vec_complex[2] = std::complex<double>(a.tau[j][0][2], 0.);

    K_rot_rot(j_vec_complex, a.rkt[i], a.root[j], a.norm[j], integral_par_f_simple,
            integral_par_f_grad_simple, param, param_seg, cur_res3);
    block2[0][0] = scal_prod(cur_res3, a.tau[i][0]);
    block2[1][0] = scal_prod(cur_res3, a.tau[i][1]);


    j_vec_complex[0] = std::complex<double>(a.tau[j][1][0], 0.);
    j_vec_complex[1] = std::complex<double>(a.tau[j][1][1], 0.);
    j_vec_complex[2] = std::complex<double>(a.tau[j][1][2], 0.);

    K_rot_rot(j_vec_complex, a.rkt[i], a.root[j], a.norm[j], integral_par_f_simple, integral_par_f_grad_simple, param, param_seg, cur_res3);
    block2[0][1] = scal_prod(cur_res3, a.tau[i][0]);
    block2[1][1] = scal_prod(cur_res3, a.tau[i][1]);

    return;

}


void get_block2_R(const TGrid_DC_Full& a, int i, int j,
        integral_par int_par_f_R,
        f_par param_f_grad_R, std::complex<double> (*block2)[2],
        int sign)
{
    if (i == j)
    {
//R = 0
        double v[3];
        block2[0][0] = 0.;
        block2[1][1] = 0.;

        vec_prod(a.tau[i][0], a.norm[i], v); //m = 0
        block2[1][0] = sign * 0.5 * scal_prod(v, a.tau[i][1]);

        vec_prod(a.tau[i][1], a.norm[i], v); //m = 1
        block2[0][1] = sign * 0.5 * scal_prod(v, a.tau[i][0]);

    } else {

        std::complex<double> cur_res3[3];
        std::complex<double> j_vec_complex[3];


        j_vec_complex[0] = std::complex<double>(a.tau[j][0][0], 0.);
        j_vec_complex[1] = std::complex<double>(a.tau[j][0][1], 0.);
        j_vec_complex[2] = std::complex<double>(a.tau[j][0][2], 0.);

        R_rot(j_vec_complex, a.rkt[i], a.root[j], int_par_f_R, param_f_grad_R, cur_res3);


        block2[0][0] = -scal_prod(cur_res3, a.tau[i][0]);
        block2[1][0] = -scal_prod(cur_res3, a.tau[i][1]);


        j_vec_complex[0] = std::complex<double>(a.tau[j][1][0], 0.);
        j_vec_complex[1] = std::complex<double>(a.tau[j][1][1], 0.);
        j_vec_complex[2] = std::complex<double>(a.tau[j][1][2], 0.);


        R_rot(j_vec_complex, a.rkt[i], a.root[j], int_par_f_R, param_f_grad_R, cur_res3);

        block2[0][1] = -scal_prod(cur_res3, a.tau[i][0]);
        block2[1][1] = -scal_prod(cur_res3, a.tau[i][1]);
    }
    return;
}

