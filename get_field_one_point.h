#ifndef _GET_FIELD_ONE_POINT_H_
#define _GET_FILED_ONE_POINT_H_

#include <complex>
#include "common_type.h"
#include "ED_Par.h"
#include "Num_Par.h"

//===================================================Подсчет токов==================================
void get_field_area1(std::complex<double>** j_E, std::complex<double>** j_M,
        const TGrid_DC_Full& a,
        ED_Par ed_param, Num_Par num_param,
        double* e0,
        double* x, std::complex<double>* field_val);

void get_field_area2(std::complex<double>** j_E, std::complex<double>** j_M,
        const TGrid_DC_Full& a,
        ED_Par ed_param, Num_Par num_param,
        double* e0,
        double* x, std::complex<double>* field_val);
#endif

