#include <iostream>
#include <cmath>
#include <complex>

//===========================================Absolute value=========================================

double abs_tmp(const double a)
{
    return fabs(a);
}

double abs_tmp(const std::complex<double> a)
{
    return sqrt(a.real() * a.real() + a.imag() * a.imag());
}

//===========================================My comparate===========================================

int cmp(const std::complex<double> a, const std::complex<double> b)
{
    double a_abs = sqrt(a.real() * a.real() + a.imag() * a.imag());
    double b_abs = sqrt(b.real() * b.real() + b.imag() * b.imag());
    if (a_abs < b_abs)
    {
        return -1;
    } else if (a_abs > b_abs)
    {
        return 1;
    } else return 0;
}

int cmp(const std::complex<double> a, const double b)
{
    double a_abs = sqrt(a.real() * a.real() + a.imag() * a.imag());
    if (a_abs < b)
    {
        return -1;
    } else if (a_abs > b)
    {
        return 1;
    } else return 0;
}

int cmp(const double a, const std::complex<double> b)
{
    double b_abs = sqrt(b.real() * b.real() + b.imag() * b.imag());
    if (a < b_abs)
    {
        return -1;
    } else if (a > b_abs)
    {
        return 1;
    } else return 0;
}

int cmp(const double a, const double b)
{
    if (a < b)
    {
        return -1;
    } else if (a > b)
    {
        return 1;
    } else return 0;
}


//=======================================Scalar product=============================================

double scal_prod(const double* vec_1, const double* vec_2)
{
    return (vec_1[0] * vec_2[0] + vec_1[1] * vec_2[1] + vec_1[2] * vec_2[2]);
}

std::complex<double> scal_prod(const std::complex<double>* vec_1, const double* vec_2)
{
    return (vec_1[0] * vec_2[0] + vec_1[1] * vec_2[1] + vec_1[2] * vec_2[2]);
}

std::complex<double> scal_prod(const double* vec_1, const std::complex<double>* vec_2)
{
    return (vec_1[0] * vec_2[0] + vec_1[1] * vec_2[1] + vec_1[2] * vec_2[2]);
}

std::complex<double> scal_prod(const std::complex<double>* vec_1, const std::complex<double>* vec_2)
{
    return (vec_1[0] * vec_2[0] + vec_1[1] * vec_2[1] + vec_1[2] * vec_2[2]);
}

//=======================================Vector product=============================================
void vec_prod(const double* vec_1, const double* vec_2, double* res)
{
    res[0] = vec_1[1] * vec_2[2] - vec_1[2] * vec_2[1];
    res[1] = vec_1[2] * vec_2[0] - vec_1[0] * vec_2[2];
    res[2] = vec_1[0] * vec_2[1] - vec_1[1] * vec_2[0];
    return;
}

void vec_prod(const std::complex<double>* vec_1, const double* vec_2, std::complex<double>* res)
{
    res[0] = vec_1[1] * vec_2[2] - vec_1[2] * vec_2[1];
    res[1] = vec_1[2] * vec_2[0] - vec_1[0] * vec_2[2];
    res[2] = vec_1[0] * vec_2[1] - vec_1[1] * vec_2[0];
    return;
}

void vec_prod(const double* vec_1, const std::complex<double>* vec_2, std::complex<double>* res)
{
    res[0] = vec_1[1] * vec_2[2] - vec_1[2] * vec_2[1];
    res[1] = vec_1[2] * vec_2[0] - vec_1[0] * vec_2[2];
    res[2] = vec_1[0] * vec_2[1] - vec_1[1] * vec_2[0];
    return;
}

void vec_prod(const std::complex<double>* vec_1, const std::complex<double>* vec_2, std::complex<double>* res)
{
    res[0] = vec_1[1] * vec_2[2] - vec_1[2] * vec_2[1];
    res[1] = vec_1[2] * vec_2[0] - vec_1[0] * vec_2[2];
    res[2] = vec_1[0] * vec_2[1] - vec_1[1] * vec_2[0];
    return;
}

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
double solid_angle(const double* x_a, const double* x_b, const double* x_c, const double* x)
{
    double r1[3], r2[3], r3[3], l1, l2, l3, deg, v[3];
    r1[0] = x_a[0] - x[0];
    r1[1] = x_a[1] - x[1];
    r1[2] = x_a[2] - x[2];

    r2[0] = x_b[0] - x[0];
    r2[1] = x_b[1] - x[1];
    r2[2] = x_b[2] - x[2];

    r3[0] = x_c[0] - x[0];
    r3[1] = x_c[1] - x[1];
    r3[2] = x_c[2] - x[2];

    l1 = vec_length(r1);
    l2 = vec_length(r2);
    l3 = vec_length(r3);

    vec_prod(r2, r3, v);

    deg = scal_prod(r1, v) / (l1 * l2 * l3 +
            scal_prod(r1, r2) * l3 + scal_prod(r2, r3) * l1 + scal_prod(r3, r1) * l2);
    return 2 * atan(deg);
}

//==========================================Cell normal=============================================
void norm_func(double** rut0, double* norm_res)
{
    double ac[3], bd[3];
    for (int i = 0; i < 3; i++)
    {
        ac[i] = rut0[2][i] - rut0[0][i];
        bd[i] = rut0[3][i] - rut0[1][i]; 
    }
    double v[3], l;
    vec_prod(ac, bd, v);
    l = vec_length(v);
    norm_res[0] = v[0] / l;
    norm_res[1] = v[1] / l;
    norm_res[2] = v[2] / l;
}
