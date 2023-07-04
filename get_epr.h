#ifndef _GET_EPR_H_
#define _GET_EPR_H_

#include <complex>
#include <string>
#include "ED_Par.h"
#include "common_type.h"

int get_epr_ideal_conductor(double* epr, double* epr_ln,
        std::complex<double>** j,
        ED_Par ed_param, const TGrid_DC_Full& a, double* e0, std::string filename);

int get_epr_dielectric(double* epr, double* epr_ln,
        std::complex<double>** j_E, std::complex<double>** j_M,
        ED_Par ed_param, const TGrid_DC_Full& a, double* e0, std::string filename);

#endif
