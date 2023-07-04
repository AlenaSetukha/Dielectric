#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>
#include <string>

#include "element_geom.h"
#include "common_type.h"
#include "ED_Par.h"
#include "Num_Par.h"
#include "constants.h"
#include "e0.h"
#include "get_area.h"
#include "get_field_one_point.h"


//filename_in - имя файла выхода(без суффиксов)
//a - класс, описывающий сетку
//j_E - вектор с электрическими токами
//j_M - вектор с магнитными токами
//ed_param - э/д параметры(каждой области)
//num_param - численные параметры для инетгралов
//polar_type - 0 - H(горизонтальная поляризация), 1 - V(вертикальная поляризация)
//filename - файл с точками, в которых считаем поле



int get_field(std::string filename_in,
        const TGrid_DC_Full& a,
        std::complex<double>** j_E , std::complex<double>** j_M,
        ED_Par ed_param, Num_Par num_param,
        int polar_type, std::string filename)
{
//=============================Создание файлов===================================
    std::ofstream fout_u_real(filename_in + "_real.gv");
    if (!fout_u_real.is_open()) {
        std::cout << "Open " << filename_in << "_real.gv error" << std::endl;
        return -1;
    }

    std::ofstream fout_u_image(filename_in + "_image.gv");
    if (!fout_u_image.is_open()) {
        std::cout << "Open " << filename_in << "_image.gv error" << std::endl;
        return -1;
    }

    std::ofstream fout_u_abs(filename_in + "_abs.gdr"); //Модуль электрическое поля
    if (!fout_u_abs.is_open()) {
        std::cout << "Open " << filename_in << "_abs.gdr error" << std::endl;
        return -1;
    }

    std::ifstream fin(filename);                // Точки, в которых считаем поле
    if (!fin.is_open()) {
        std::cout << "Open " << filename << " error" << std::endl;
        return -1;
    }

//========================Расчет поля в точках из файла===========================

    int n, n1, n2, area;
    fin >> n1 >> n2;


    fout_u_real << n1 << " " << n2 << std::endl;
    fout_u_image << n1 << " " << n2 << std::endl;
    fout_u_abs << n1 << " " << n2 << std::endl;

    E0_Polar_Type e0(polar_type, ed_param.k_vec);
    Constants c;
    n = n1 * n2;//количество точек, в которых хотим считать поле
    double x[3]; // точка, в которой считаем поле
    std::complex<double> field_val[3];

    for (int i = 0; i < n; i++)
    {
std::cout << i << std::endl;
        fin >> x[0] >> x[1] >> x[2];
//Определение, в какой области лежит точка(вне или внутри диэлектрика)
        area = get_area(x, a);
        if (area == -1)
        {
            get_field_area2(j_E, j_M, a, ed_param, num_param, e0.e0, x, field_val);
        } else {
            get_field_area1(j_E, j_M, a, ed_param, num_param, e0.e0, x, field_val);
        }
        fout_u_real << field_val[0].real() << " " << field_val[1].real() << " " << field_val[2].real() << std::endl;
        fout_u_image << field_val[0].imag() << " " << field_val[1].imag() << " " << field_val[2].imag() << std::endl;
        fout_u_abs << vec_length(field_val) << std::endl;
    }



//Закрытие файлов
    fout_u_real.close();
    fout_u_image.close();
    fout_u_abs.close();
    fin.close();
    return 0;
}








int get_field_exact(std::string filename_in,
        const TGrid_DC_Full& a,
        ED_Par ed_param,
        int polar_type, std::string filename)
{
//=============================Создание файлов===================================
    std::ofstream fout_u_real(filename_in + "_real.gv");
    if (!fout_u_real.is_open()) {
        std::cout << "Open " << filename_in << "_real.gv error" << std::endl;
        return -1;
    }

    std::ofstream fout_u_image(filename_in + "_image.gv");
    if (!fout_u_image.is_open()) {
        std::cout << "Open " << filename_in << "_image.gv error" << std::endl;
        return -1;
    }

    std::ofstream fout_u_abs(filename_in + "_abs.gdr"); //Модуль электрическое поля
    if (!fout_u_abs.is_open()) {
        std::cout << "Open " << filename_in << "_abs.gdr error" << std::endl;
        return -1;
    }

    std::ifstream fin(filename);                // Точки, в которых считаем поле
    if (!fin.is_open()) {
        std::cout << "Open " << filename << " error" << std::endl;
        return -1;
    }

//========================Расчет поля в точках из файла===========================

    int n, n1, n2, area;
    fin >> n1 >> n2;


    fout_u_real << n1 << " " << n2 << std::endl;
    fout_u_image << n1 << " " << n2 << std::endl;
    fout_u_abs << n1 << " " << n2 << std::endl;

    E0_Polar_Type e0(polar_type, ed_param.k_vec);
    Constants c;
    n = n1 * n2;//количество точек, в которых хотим считать поле
    double x[3]; // точка, в которой считаем поле
    std::complex<double> field_val[3], deg;

    for (int i = 0; i < n; i++)
    {
        fin >> x[0] >> x[1] >> x[2];
//Определение, в какой области лежит точка
        area = get_area(x, a);
        if (area == 1)//поле снаружи = 0
        {
            field_val[0] = std::complex<double>(0., 0.);
            field_val[1] = std::complex<double>(0., 0.);
            field_val[2] = std::complex<double>(0., 0.);
        } else {
            deg = pow(c.e, c.i_complex * scal_prod(ed_param.k_vec, x));
            field_val[0] = e0.e0[0] * deg;
            field_val[1] = e0.e0[1] * deg;
            field_val[2] = e0.e0[2] * deg;
        }
        fout_u_real << field_val[0].real() << " " << field_val[1].real() << " " << field_val[2].real() << std::endl;
        fout_u_image << field_val[0].imag() << " " << field_val[1].imag() << " " << field_val[2].imag() << std::endl;
        fout_u_abs << vec_length(field_val) << std::endl;
    }



//Закрытие файлов
    fout_u_real.close();
    fout_u_image.close();
    fout_u_abs.close();
    fin.close();
    return 0;
}







int get_field_omega2(std::string filename_in,
        const TGrid_DC_Full& a,
        std::complex<double>** j_E , std::complex<double>** j_M,
        ED_Par ed_param, Num_Par num_param,
        int polar_type, std::string filename)
{
//=============================Создание файлов===================================
    std::ofstream fout_u_real(filename_in + "_real.gv");
    if (!fout_u_real.is_open()) {
        std::cout << "Open " << filename_in << "_real.gv error" << std::endl;
        return -1;
    }

    std::ofstream fout_u_image(filename_in + "_image.gv");
    if (!fout_u_image.is_open()) {
        std::cout << "Open " << filename_in << "_image.gv error" << std::endl;
        return -1;
    }

    std::ofstream fout_u_abs(filename_in + "_abs.gdr"); //Модуль электрическое поля
    if (!fout_u_abs.is_open()) {
        std::cout << "Open " << filename_in << "_abs.gdr error" << std::endl;
        return -1;
    }

    std::ifstream fin(filename);                // Точки, в которых считаем поле
    if (!fin.is_open()) {
        std::cout << "Open " << filename << " error" << std::endl;
        return -1;
    }

//========================Расчет поля в точках из файла===========================

    int n, n1, n2;
    fin >> n1 >> n2;


    fout_u_real << n1 << " " << n2 << std::endl;
    fout_u_image << n1 << " " << n2 << std::endl;
    fout_u_abs << n1 << " " << n2 << std::endl;

    E0_Polar_Type e0(polar_type, ed_param.k_vec);
    Constants c;
    n = n1 * n2;//количество точек, в которых хотим считать поле
    double x[3]; // точка, в которой считаем поле
    std::complex<double> field_val[3];

    for (int i = 0; i < n; i++)
    {
        fin >> x[0] >> x[1] >> x[2];
//Определение, в какой области лежит точка
        field_val[0] = std::complex<double>(0., 0.);
        field_val[1] = std::complex<double>(0., 0.);
        field_val[2] = std::complex<double>(0., 0.);
        get_field_area2(j_E, j_M, a, ed_param, num_param, e0.e0, x, field_val);
        fout_u_real << field_val[0].real() << " " << field_val[1].real() << " " << field_val[2].real() << std::endl;
        fout_u_image << field_val[0].imag() << " " << field_val[1].imag() << " " << field_val[2].imag() << std::endl;
        fout_u_abs << vec_length(field_val) << std::endl;
    }



//Закрытие файлов
    fout_u_real.close();
    fout_u_image.close();
    fout_u_abs.close();
    fin.close();
    return 0;
}



