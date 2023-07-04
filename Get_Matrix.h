#ifndef _GET_MATRIX_H_
#define _GET_MATRIX_H_

#include <complex>
#include "common_type.h"
#include "ED_Par.h"
#include "Num_Par.h"

//=================================================Набор матрицы A=====================================================
void get_matrix_1block(std::complex<double>** matrix, const TGrid_DC_Full& a, ED_Par ed_param, Num_Par num_param);


//=================================================Набор матрицы A|B|C|D=====================================================
void get_matrix_4block(std::complex<double>** matrix, const TGrid_DC_Full& a, ED_Par ed_param, Num_Par num_param);


#endif
