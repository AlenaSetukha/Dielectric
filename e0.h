#ifndef _POLAR_TYPE_H_
#define _POLAR_TYPE_H_

class E0_Polar_Type
{
public:
    int type;//тип поляризации
    double e0[3];//вектор поляризации

    E0_Polar_Type(int type_in, double* k_in);

};
#endif

