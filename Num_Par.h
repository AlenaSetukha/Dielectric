#ifndef _NUM_PAR_H_
#define _NUM_PAR_H_

class Num_Par
{
public:
    double eps, rs;// точность расчета интегралов, радиус сглаживания
    int n_start_seg, n_start, p_max, p_max_seg;//стартовое разбиение на отрезке, на ячейке, предельное разбиение на ячейке, отрезке

    Num_Par(const char* filename);
};
#endif
