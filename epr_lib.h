#ifndef _EPR_LIB_H_
#define _EPR_LIB_H_

#include <complex>
#include "common_type.h"
#include "ED_Par.h"

//------------------------------------Функцияя эпр для идеального проводника-------------------------------------
double f_epr_ideal_conductor(double* tau, const TGrid_DC_Full& a,
        std::complex<double> k, std::complex<double>** j, double* e0);

//------------------------------------Функцияя эпр для диэлектрика-------------------------------------
double f_epr_dielectric(double* tau, const TGrid_DC_Full& a,
        std::complex<double>** j_E, std::complex<double>** j_M,
        double* e0, ED_Par ed_param);
#endif