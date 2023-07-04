#ifndef _ELEMENT_GEOM_H_
#define _ELEMENT_GEOM_H_

#include <complex>

//===========================================Absolute value=========================================

double abs_tmp(const double a);

double abs_tmp(const std::complex<double> a);

//===========================================My comparate===========================================

int cmp(const std::complex<double> a, const std::complex<double> b);

int cmp(const std::complex<double> a, const double b);

int cmp(const double a, const std::complex<double> b);

int cmp(const double a, const double b);


//=======================================Scalar product=============================================

double scal_prod(const double* vec_1, const double* vec_2);

std::complex<double> scal_prod(const std::complex<double>* vec_1, const double* vec_2);

std::complex<double> scal_prod(const double* vec_1, const std::complex<double>* vec_2);

std::complex<double> scal_prod(const std::complex<double>* vec_1, const std::complex<double>* vec_2);

//=======================================Vector product=============================================
void vec_prod(const double* vec_1, const double* vec_2, double* res);

void vec_prod(const std::complex<double>* vec_1, const double* vec_2, std::complex<double>* res);

void vec_prod(const double* vec_1, const std::complex<double>* vec_2, std::complex<double>* res);

void vec_prod(const std::complex<double>* vec_1, const std::complex<double>* vec_2, std::complex<double>* res);


//=====================================Distance btw 2 points========================================
template <typename T>
T dist(const T* vec_1, const T* vec_2)
{
    return sqrt((vec_1[0] - vec_2[0]) * (vec_1[0] - vec_2[0]) +
        (vec_1[1] - vec_2[1]) * (vec_1[1] - vec_2[1]) +
        (vec_1[2] - vec_2[2]) * (vec_1[2] - vec_2[2]));
}

//=======================================Square of triagnle=========================================
template <typename T>
T tr_square(const T* pnt_1, const T* pnt_2, const T* pnt_3)
{
    T ab[3], ac[3];
    for (int i = 0; i < 3; i++) {
        ab[i] = pnt_2[i] - pnt_1[i];
        ac[i] = pnt_3[i] - pnt_1[i];
    }
    T* n = vec_prod(ab, ac);
    T res = scal_prod(n, n);
    res = sqrt(res) / (static_cast<T>(2));
    return res;
}

//=========================================Vector length============================================
template <typename T>
double vec_length(const T* vec_1)
{
    return sqrt(abs_tmp(vec_1[0]) * abs_tmp(vec_1[0]) + abs_tmp(vec_1[1]) *
                    abs_tmp(vec_1[1]) + abs_tmp(vec_1[2]) * abs_tmp(vec_1[2]));
}

//==========================================Solid angle=============================================
double solid_angle(const double* x_a, const double* x_b, const double* x_c, const double* x);

//==========================================Cell normal=============================================
void norm_func(double** rut0, double* norm_res);

#endif
