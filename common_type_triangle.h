#ifndef _COMMON_TYPE_TRIANGLE_H_
#define _COMMON_TYPE_TRIANGLE_H_

#include <string>

class TGrid_DC_Full_Triangle
{
public:
    int num_obj;
    int num_mod;
    int num_frm;

    int** beg_end;// [mod][2] номер 1 и посл треугольника на модуле с номером mod
    double*** root;// [frm][3][3]
    double* s;//[frm] площадь ячейки
    double** norm;// [frm][3]
    double** rkt;// [frm][3] - точки коллокации(центр ячейки)
    double*** tau;// [frm][2][3] - касательные векторы к ячейке
    int** beg_end_obj;//[obj][2] 1 и посл модуль на каждом объекте с нoмером obj(сплошная сетка)

    double max_diag;


    void fill_TGrid_Triangle(std::string filename);
    //------------------------------------------Constructor-------------------------------------------------------
    TGrid_DC_Full_Triangle(std::string filename);//сквозная сетка на все объекты и все модули вместе
    //----------------------------------------Copy Constructor----------------------------------------------------
    TGrid_DC_Full_Triangle(const TGrid_DC_Full_Triangle &obj);
    //------------------------------------------Destructor--------------------------------------------------------
    ~TGrid_DC_Full_Triangle();

};
#endif
