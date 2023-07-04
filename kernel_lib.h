#ifndef _KERNEL_LIB_H_
#define _KERNEL_LIB_H_

#include <complex>
#include "f_par.h"

//-----------------------Потенциал простого слоя для ур-я Лапласа----------------------
void f_simple_pot_L(double* x, double* y, f_par param, double* ff);

//---------------------Потенциал простого слоя для ур-я Гельмгольца--------------------
void f_simple_pot_G(double* x, double* y, f_par param, std::complex<double>* ff);

//-----------------------Потенциал двойного слоя для ур-я Гельмгольца----------------------
void f_double_pot_G(double* x, double* y, f_par param,  std::complex<double>* ff);

//---------------------Потенциал двойного слоя для ур-я Лапласа--------------------+++
void f_double_pot_L(double* x, double* y, f_par param,  double* ff);

//---------------------Векторный потенциал для уравнения Лапласа-------------------
template <typename P>
void f_vector_pot_L(double* x, double* y, f_par param, P* ff);

//---------------Градиент потенциала простого слоя уравнения Гельмгольца---------------
void f_grad_simple_pot_G(double* x, double* y, f_par param, std::complex<double>* ff);

//---------------Градиент потенциала простого слоя уравнения Гельмгольца eps-----------
void f_grad_simple_pot_G_eps(double* x, double* y, f_par param, std::complex<double>* ff);

#endif

