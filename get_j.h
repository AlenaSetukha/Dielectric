#ifndef _GET_J_H_
#define _GET_J_H_

#include <complex>
#include <string>
#include "common_type.h"

//===================================================Подсчет токов==================================
int  get_j(std::complex<double>**  j_vec, const TGrid_DC_Full& a, std::complex<double>* b, std::string filename);

int get_j_from_files(std::complex<double>**  j_vec, std::string filename_real, std::string filename_image);
#endif

