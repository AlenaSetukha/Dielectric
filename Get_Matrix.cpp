#include <iostream>
#include <cmath>
#include <complex>

#include "common_type.h"
#include "ED_Par.h"
#include "Num_Par.h"
#include "integral_par.h"
#include "f_par.h"
#include "element_geom.h"
#include "block.h"
#include "constants.h"

//=================================================Набор матрицы A=====================================================
void get_matrix_1block(std::complex<double>** matrix, const TGrid_DC_Full& a, ED_Par ed_param, Num_Par num_param)
{
//========================Инициализация параметров===========================

//-----------f_simple_pot_G---
    integral_par f_simple_pot_G_par(1, num_param.n_start, num_param.p_max, num_param.eps);//1,10,1,0.001
    f_par param(a.max_diag * num_param.rs, ed_param.k[0]);

//-----------f_grad_simple_pot_G---
    integral_par f_grad_simple_pot_G_par(3, num_param.n_start_seg, num_param.p_max_seg, num_param.eps);//3,10,1,0.001
    f_par param_seg(0.000000000001, ed_param.k[0]);

//==========================Создание СЛАУ====================
    int i0, j0;
    std::complex<double> block[2][2];

    for (int i = 0; i < a.num_frm; i++)
    {
        std::cout << " i " << i << std::endl;
        i0 = 2 * i;
        for (int j = 0; j < a.num_frm; j++)
        {
            j0 = 2 * j;
            get_block2_K(a, i, j, f_simple_pot_G_par, f_grad_simple_pot_G_par, param, param_seg, block);
            matrix[i0][j0] = block[0][0];
            matrix[i0][j0 + 1] = block[0][1];
            matrix[i0 + 1][j0] = block[1][0];
            matrix[i0 + 1][j0 + 1] = block[1][1];
        }
    }
    return;
}





//=================================================Набор матрицы A|B|C|D=====================================================
void get_matrix_4block(std::complex<double>** matrix, const TGrid_DC_Full& a, ED_Par ed_param, Num_Par num_param)
{
    Constants c;
//==========================Заполнение СЛАУ====================
    int i0, j0;
    std::complex<double> block[2][2];

//========================Инициализация параметров=============
    integral_par int_par_f_grad_K(3, num_param.n_start_seg, num_param.p_max_seg, num_param.eps);//размернсть, стартовое разбиение, максимальное разбиение, точность
    integral_par int_par_f_grad_R(3, num_param.n_start, num_param.p_max, num_param.eps);
    integral_par int_par_f_K(1, num_param.n_start, num_param.p_max, num_param.eps);
    
    
    f_par param_f_grad_K(0.000000000001, ed_param.k[0]);
    f_par param_f_K(a.max_diag * num_param.rs, ed_param.k[0]);
    f_par param_f_grad_R(a.max_diag * num_param.rs, ed_param.k[0]);

    std::complex<double> mult1 = c.i_complex / ed_param.omega / ed_param.eps_d[0] / c.eps0;//множитель
    std::complex<double> mult2 = c.i_complex / ed_param.omega / ed_param.eps_d[1] / c.eps0;

//=======================Набор блока A=========================
    for (int i = 0; i < a.num_frm; i++)
    {
        std::cout << " i " << i << std::endl;
        i0 = 2 * i;
        for (int j = 0; j < a.num_frm; j++)
        {
            j0 = 2 * j;
            get_block2_K(a, i, j, int_par_f_K, int_par_f_grad_K, param_f_K, param_f_grad_K, block);
            matrix[i0][j0] = block[0][0] * mult1;
            matrix[i0][j0 + 1] = block[0][1] * mult1;
            matrix[i0 + 1][j0] = block[1][0] * mult1;
            matrix[i0 + 1][j0 + 1] = block[1][1] * mult1;
        }
    }
//=======================Набор блока B=========================

    for (int i = 0; i < a.num_frm; i++)
    {
        std::cout << " i " << i << std::endl;
        i0 = 2 * i;
        for (int j = 0; j < a.num_frm; j++)
        {
            j0 = 2 * j + 2 * a.num_frm;
            get_block2_R(a, i, j, int_par_f_grad_R, param_f_grad_R, block, 1);
            matrix[i0][j0] = block[0][0];
            matrix[i0][j0 + 1] = block[0][1];
            matrix[i0 + 1][j0] = block[1][0];
            matrix[i0 + 1][j0 + 1] = block[1][1];
        }
    }



//=======================Набор блока C=========================
    param_f_grad_K.k = ed_param.k[1];
    param_f_K.k = ed_param.k[1];

    for (int i = 0; i < a.num_frm; i++)
    {
        std::cout << " i " << i << std::endl;
        i0 = 2 * i + 2 * a.num_frm;
        for (int j = 0; j < a.num_frm; j++)
        {
            j0 = 2 * j;
            get_block2_K(a, i, j, int_par_f_K, int_par_f_grad_K, param_f_K, param_f_grad_K, block);
            matrix[i0][j0] = block[0][0] * mult2;
            matrix[i0][j0 + 1] = block[0][1] * mult2;
            matrix[i0 + 1][j0] = block[1][0] * mult2;
            matrix[i0 + 1][j0 + 1] = block[1][1] * mult2;
        }
    }




//=======================Набор блока D=========================
    param_f_grad_R.k = ed_param.k[1];

    for (int i = 0; i < a.num_frm; i++)
    {
        std::cout << " i " << i << std::endl;
        i0 = 2 * i + 2 * a.num_frm;
        for (int j = 0; j < a.num_frm; j++)
        {
            j0 = 2 * j + 2 * a.num_frm;
            get_block2_R(a, i, j, int_par_f_grad_R, param_f_grad_R, block, -1);
            matrix[i0][j0] = block[0][0];
            matrix[i0][j0 + 1] = block[0][1];
            matrix[i0 + 1][j0] = block[1][0];
            matrix[i0 + 1][j0 + 1] = block[1][1];
        }
    }
    return;
}
