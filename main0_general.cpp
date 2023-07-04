#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>


#include "constants.h"
#include "common_type.h"
#include "ED_Par.h"
#include "Num_Par.h"
#include "Get_Matrix.h"
#include "e0.h"
#include "get_b.h"
#include "get_j.h"
#include "get_field.h"
#include "get_epr.h"
#include "solve.h"

int main()
{
//======================Constants==============================
    Constants c;
//======================Геометрия==============================
    TGrid_DC_Full a("data//geodat.dat");
    a.fill_TGrid("data//geodat.dat");

//======================Parameters=============================
    ED_Par ed_param("data//ed_param.txt");
    Num_Par num_param("data//num_param.txt");

//E0
    E0_Polar_Type e0_H(0, ed_param.k_vec);
    E0_Polar_Type e0_V(1, ed_param.k_vec);

//======================Запись параметров в файл=============================
    std::ofstream fout_params("my_results//ed_params_res.txt");
    if (!fout_params.is_open()) {
        std::cout << "Open my_results//ed_params_res.txt error" << std::endl;
        return -1;
    }

    fout_params << "k_inc = " <<  ed_param.k[0] << std::endl;
    fout_params << "omega = " << ed_param.omega << std::endl;
    fout_params << "num_frm = " <<  a.num_frm << std::endl;
    fout_params << "partition diameter = " << a.max_diag << std::endl;
    fout_params << "k_1 = " <<  ed_param.k[0] << std::endl;
    fout_params<< "k_2 = " <<  ed_param.k[1] << std::endl;
    fout_params << "eps1 = " <<  ed_param.eps_d[0] * c.eps0 << std::endl;
    fout_params << "eps2 = " <<  ed_param.eps_d[1] * c.eps0<< std::endl;
    fout_params << "mu1 = " <<  ed_param.m_d[0] * c.m0 << std::endl;
    fout_params << "mu2 = " <<  ed_param.m_d[1] * c.m0 << std::endl;
    fout_params << "lambda1 = " <<  ed_param.lambda[0] << std::endl;
    fout_params << "lambda2 = " <<  ed_param.lambda[1] << std::endl;
    fout_params << "k_vec: " << ed_param.k_vec[0] << " " << ed_param.k_vec[1] << " " << ed_param.k_vec[2] << std::endl;
    fout_params << "e0_H: " << e0_H.e0[0] << " " << e0_H.e0[1] << " " << e0_H.e0[2] << std::endl;
    fout_params << "e0_V: " << e0_V.e0[0] << " " << e0_V.e0[1] << " " << e0_V.e0[2] << std::endl;

    fout_params.close();
//=====================Создание СЛАУ===========================
    int N = a.num_frm * 4;
    std::complex<double>** matrix = new std::complex<double>*[N];
    for (int i = 0; i < N; i++)
    {
        matrix[i] = new std::complex<double>[N];
    }

    get_matrix_4block(matrix, a, ed_param, num_param);
std::cout << "matrix4block done" << std::endl;
//===================Создание правой части=====================
    std::complex<double>* b_H = new std::complex<double>[N];//4n
    std::complex<double>* b_V = new std::complex<double>[N];//4n

    get_b(b_H, e0_H.e0, ed_param.k_vec, a);//заполнит 2n, остальное нули
    get_b(b_V, e0_V.e0, ed_param.k_vec, a);//заполнит 2n, остальное нули

    for (int i = N/2; i < N; i++)
    {
        b_H[i] = std::complex<double>(0., 0.);
        b_V[i] = std::complex<double>(0., 0.);
    }
//========================Решение СЛАУ==========================
    int* ip = new int[N];
    decomp(N, matrix, ip);

    solve(N, matrix, b_H, ip);// b_H - поле в точках коллокации при горизонтальной поляризации; первые 2n чисел - je, оставшиеся 2n - jm
    solve(N, matrix, b_V, ip);// b_V - поле в точках коллокации при вертикальной поляризации; первые 2n чисел - je, оставшиеся 2n - jm
std::cout << "solve done!" << std::endl;
//=====================Запись точек коллокации в файл===========
    std::ofstream fout("my_results//body.gr");
    if (!fout.is_open()) {
        std::cout << "Open my_results//body.gr error" << std::endl;
        return -1;
    }
    fout << 2 << " " << a.num_frm / 2 << std::endl;

    for (int i = 0; i < a.num_frm; i++)
    {
        fout << a.rkt[i][0] << " " << a.rkt[i][1] << " " << a.rkt[i][2] << std::endl;
    }
    fout.close();
//=======================Нахождение токов H,V===================
    std::complex<double>** j_vec_H_E = new std::complex<double>*[a.num_frm];
    for (int i = 0; i < a.num_frm; i++)
    {
        j_vec_H_E[i] = new std::complex<double>[3];
    }

    std::complex<double>** j_vec_H_M = new std::complex<double>*[a.num_frm];
    for (int i = 0; i < a.num_frm; i++)
    {
        j_vec_H_M[i] = new std::complex<double>[3];
    }

    std::complex<double>** j_vec_V_E = new std::complex<double>*[a.num_frm];
    for (int i = 0; i < a.num_frm; i++)
    {
        j_vec_V_E[i] = new std::complex<double>[3];
    }

    std::complex<double>** j_vec_V_M = new std::complex<double>*[a.num_frm];
    for (int i = 0; i < a.num_frm; i++)
    {
        j_vec_V_M[i] = new std::complex<double>[3];
    }


    get_j(j_vec_H_E, a, b_H, "my_results//j_H_E");                         //сразу с записью в файлы "j_H_E_real.gv, j_H_E_image.gv"
    get_j(j_vec_H_M, a, &b_H[2 * a.num_frm], "my_results//j_H_M");         //сразу с записью в файлы "j_H_M_real.gv, j_H_M_image.gv"
    get_j(j_vec_V_E, a, b_V, "my_results//j_V_E");                         //сразу с записью в файлы "j_V_E_real.gv, j_V_E_image.gv"
    get_j(j_vec_V_M, a, &b_V[2 * a.num_frm], "my_results//j_V_M");         //сразу с записью в файлы "j_V_E_real.gv, j_V_E_image.gv"
std::cout << "get_j done!" << std::endl;
//==========Нахождение поля в заданных точках======================

    get_field("my_results//field_H", a, j_vec_H_E, j_vec_H_M, ed_param, num_param, 0, "data//grid.gr");//сразу с записью в файлы "field_H_real.gv", "field_H_image.gv", "field_H_abs.gdr"
    get_field("my_results//field_V", a, j_vec_V_E, j_vec_V_M, ed_param, num_param, 1, "data//grid.gr");//сразу с записью в файлы "field_V_real.gv", "field_V_image.gv", "field_V_abs.gdr"

std::cout << "get_field done!" << std::endl;

//========================Нахождение ЭПР H,V=======================

    int epr_n = ed_param.alpha_end - ed_param.alpha_start;//180

    double* epr_H = new double[epr_n];
    double* epr_V = new double[epr_n];
    double* epr_ln_H = new double[epr_n];
    double* epr_ln_V = new double[epr_n];
    get_epr_dielectric(epr_H, epr_ln_H, j_vec_H_E, j_vec_H_M, ed_param, a, e0_H.e0, "my_results//epr_h.txt");//сразу с записью в файл "epr_h.txt"
    get_epr_dielectric(epr_V, epr_ln_V, j_vec_V_E, j_vec_V_M, ed_param, a, e0_V.e0, "my_results//epr_v.txt");//сразу с записью в файл "epr_v.txt"

std::cout << "get_epr done!" << std::endl;

//=========================Очистка памяти===================
std::cout << "+++" << std::endl;

    for (int i = 0; i < N; i++) {
        delete[] matrix[i];
    }
    delete[] matrix;

    delete[] b_H;
    delete[] b_V;
    delete[] ip;

    for (int i = 0; i < a.num_frm; i++) {
        delete[] j_vec_H_E[i];
    }
    delete[] j_vec_H_E;

    for (int i = 0; i < a.num_frm; i++) {
        delete[] j_vec_H_M[i];
    }
    delete[] j_vec_H_M;

    for (int i = 0; i < a.num_frm; i++) {
        delete[] j_vec_V_E[i];
    }
    delete[] j_vec_V_E;

    for (int i = 0; i < a.num_frm; i++) {
        delete[] j_vec_V_M[i];
    }
    delete[] j_vec_V_M;


    delete[] epr_H;
    delete[] epr_V;
    delete[] epr_ln_H;
    delete[] epr_ln_V;

    return 0;
}
