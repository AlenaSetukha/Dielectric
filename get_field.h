#ifndef _GET_FIELD_H_
#define _GET_FIELD_H_

#include <complex>
#include <string>
#include "common_type.h"

int get_field(std::string filename_in,
        const TGrid_DC_Full& a,
        std::complex<double>** j_E , std::complex<double>** j_M,
        ED_Par ed_param, Num_Par num_param,
        int polar_type, std::string filename);
int get_field_exact(std::string filename_in,
        const TGrid_DC_Full& a,
        ED_Par ed_param,
        int polar_type, std::string filename);

int get_field_omega2(std::string filename_in,
        const TGrid_DC_Full& a,
        std::complex<double>** j_E , std::complex<double>** j_M,
        ED_Par ed_param, Num_Par num_param,
        int polar_type, std::string filename);

#endif
