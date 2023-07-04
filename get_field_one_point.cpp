#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>

#include "common_type.h"
#include "ED_Par.h"
#include "Num_Par.h"
#include "integral_par.h"
#include "f_par.h"
#include "constants.h"
#include "K.h"
#include "R.h"
#include "element_geom.h"
#include "get_einc.h"


//==================================Description===============================
// Поле вычисляется в заданной точке
// В зависимости от ее положения, вызывается либо функция get_field_area1(внешняя область диэлектрика),
// либо функция get_field_aea2(внутренняя область диэлектрика).



//========================================Получение поля во внешней области==================================
void get_field_area1(std::complex<double>** j_E, std::complex<double>** j_M,
        const TGrid_DC_Full& a,
        ED_Par ed_param, Num_Par num_param,
        double* e0,
        double* x, std::complex<double>* field_val)
{
//========================Инициализация параметров===========================
    integral_par int_par_f_grad_K(3, num_param.n_start_seg, num_param.p_max_seg, num_param.eps);
    integral_par int_par_f_grad_R(3, num_param.n_start, num_param.p_max, num_param.eps);
    integral_par int_par_f_K(1, num_param.n_start, num_param.p_max, num_param.eps);

    f_par param_f_grad_K(0.000000000001, ed_param.k[0]);

    f_par param_f_K(a.max_diag * num_param.rs, ed_param.k[0]);

    f_par param_f_grad_R(a.max_diag * num_param.rs, ed_param.k[0]);

//========================Вычисление поля===========================
    Constants c;
    std::complex<double> sum_K[3], sum_R[3];
    sum_K[0] = std::complex<double>(0., 0.);
    sum_K[1] = std::complex<double>(0., 0.);
    sum_K[2] = std::complex<double>(0., 0.);
    sum_R[0] = std::complex<double>(0., 0.);
    sum_R[1] = std::complex<double>(0., 0.);
    sum_R[2] = std::complex<double>(0., 0.);
    std::complex<double> cur_res3[3];
    std::complex<double> deg = c.i_complex / (ed_param.omega * ed_param.eps_d[0] * c.eps0);
    for (int i = 0; i < a.num_frm; i++)
    {
        K_rot_rot(j_E[i], x, a.root[i], a.norm[i], int_par_f_K, int_par_f_grad_K, param_f_K, param_f_grad_K, cur_res3);
        sum_K[0] += cur_res3[0];
        sum_K[1] += cur_res3[1];
        sum_K[2] += cur_res3[2];

        R_rot(j_M[i], x, a.root[i], int_par_f_grad_R, param_f_grad_R, cur_res3);
        sum_R[0] += cur_res3[0];
        sum_R[1] += cur_res3[1];
        sum_R[2] += cur_res3[2];
    }
    std::complex<double> e_inc[3];
    get_einc(x, e0, ed_param.k_vec, e_inc);
    field_val[0] = deg * sum_K[0] - sum_R[0] + e_inc[0];
    field_val[1] = deg * sum_K[1] - sum_R[1] + e_inc[1];
    field_val[2] = deg * sum_K[2] - sum_R[2] + e_inc[2];
    return;
}





//========================================Получение поля во внутренней области==================================
void get_field_area2(std::complex<double>** j_E, std::complex<double>** j_M,
        const TGrid_DC_Full& a,
        ED_Par ed_param, Num_Par num_param,
        double* e0,
        double* x, std::complex<double>* field_val)
{
//========================Инициализация параметров===========================
    integral_par int_par_f_grad_K(3, num_param.n_start_seg, num_param.p_max_seg, num_param.eps);
    integral_par int_par_f_grad_R(3, num_param.n_start, num_param.p_max, num_param.eps);
    integral_par int_par_f_K(1, num_param.n_start, num_param.p_max, num_param.eps);

    f_par param_f_grad_K(0.000000000001, ed_param.k[1]);
    f_par param_f_K(a.max_diag * num_param.rs, ed_param.k[1]);
    f_par param_f_grad_R(a.max_diag * num_param.rs, ed_param.k[1]);
//========================Вычисление поля===========================
    Constants c;
    std::complex<double> sum_K[3], sum_R[3];
    sum_K[0] = std::complex<double>(0., 0.);
    sum_K[1] = std::complex<double>(0., 0.);
    sum_K[2] = std::complex<double>(0., 0.);
    sum_R[0] = std::complex<double>(0., 0.);
    sum_R[1] = std::complex<double>(0., 0.);
    sum_R[2] = std::complex<double>(0., 0.);
    std::complex<double> cur_res3[3];
    std::complex<double> deg = c.i_complex / (ed_param.omega * ed_param.eps_d[1] * c.eps0);

    std::complex<double> j_E_min[3], j_M_min[3];
    for (int i = 0; i < a.num_frm; i++)
    {
        j_E_min[0] = - j_E[i][0];
        j_E_min[1] = - j_E[i][1];
        j_E_min[2] = - j_E[i][2];
        K_rot_rot(j_E_min, x, a.root[i], a.norm[i], int_par_f_K, int_par_f_grad_K, param_f_K, param_f_grad_K, cur_res3);
        sum_K[0] += cur_res3[0];
        sum_K[1] += cur_res3[1];
        sum_K[2] += cur_res3[2];

        j_M_min[0] = - j_M[i][0];
        j_M_min[1] = - j_M[i][1];
        j_M_min[2] = - j_M[i][2];
        R_rot(j_M_min, x, a.root[i], int_par_f_grad_R, param_f_grad_R, cur_res3);
        sum_R[0] += cur_res3[0];
        sum_R[1] += cur_res3[1];
        sum_R[2] += cur_res3[2];
    }
    field_val[0] = deg * sum_K[0] - sum_R[0];
    field_val[1] = deg * sum_K[1] - sum_R[1];
    field_val[2] = deg * sum_K[2] - sum_R[2];
    return;
}
