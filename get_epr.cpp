#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>
#include <string>

#include "ED_Par.h"
#include "common_type.h"
#include "constants.h"
#include "epr_lib.h"


int get_epr_ideal_conductor(double* epr, double* epr_ln,
        std::complex<double>** j,
        ED_Par ed_param, const TGrid_DC_Full& a, double* e0, std::string filename)
{
//=========Создание файла для записи ЭПР====================
    std::ofstream fout_epr(filename);
    if (!fout_epr.is_open()) {
        std::cout << "Open " << filename << " error" << std::endl;
        return -1;
    }

//============Расчет ЭПР, запись в файл=====================
    Constants c;
    double tau_epr[3];

    int n = ed_param.alpha_end - ed_param.alpha_start + 1;
    for (int i = 0; i < n; i++)
    {
        tau_epr[0] = cos((ed_param.alpha_start + i) / c.ra);//перевод угла в радианы
        tau_epr[1] = sin((ed_param.alpha_start + i)  / c.ra);
        tau_epr[2] = 0.;

        epr[i] = f_epr_ideal_conductor(tau_epr, a, ed_param.k[0], j, e0);
        epr_ln[i] = 10. * log10(epr[i]);

//Угол в градусах, радианах, ЭПР, ln(ЭПР)
        fout_epr << ed_param.alpha_start + i << " " << (ed_param.alpha_start + i) / c.ra << " " << epr[i] << " " << epr_ln[i] << std::endl;
    }
//===================Закрытие файлов=======================
    fout_epr.close();
    return 0;
}




int get_epr_dielectric(double* epr, double* epr_ln,
        std::complex<double>** j_E, std::complex<double>** j_M,
        ED_Par ed_param, const TGrid_DC_Full& a, double* e0, std::string filename)
{
//=========Создание файла для записи ЭПР====================
    std::ofstream fout_epr(filename);
    if (!fout_epr.is_open()) {
        std::cout << "Open " << filename << " error" << std::endl;
        return -1;
    }

//============Расчет ЭПР, запись в файл=====================
    Constants c;
    double tau_epr[3];

    int n = ed_param.alpha_end - ed_param.alpha_start + 1;
    for (int i = 0; i < n; i++)
    {
        tau_epr[0] = cos((ed_param.alpha_start + i) / c.ra);//перевод угла в радианы
        tau_epr[1] = sin((ed_param.alpha_start + i)  / c.ra);
        tau_epr[2] = 0.;

        epr[i] = f_epr_dielectric(tau_epr, a, j_E, j_M, e0, ed_param);
        epr_ln[i] = 10. * log10(epr[i]);
//Угол в градусах, радианах, ЭПР, ln(ЭПР)
        fout_epr << ed_param.alpha_start + i << " " << (ed_param.alpha_start + i) / c.ra << " " << epr[i] << " " << epr_ln[i] << std::endl;
    }
//===================Закрытие файлов=======================
    fout_epr.close();
    return 0;
}
