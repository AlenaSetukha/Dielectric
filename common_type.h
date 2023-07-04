#ifndef _COMMON_TYPE_H_
#define _COMMON_TYPE_H_

#include <string>

class TGrid_DC_Full
{
public:
    double*** tau;// [frm][2][3] - касательные векторы к ячейке
    double** rkt;// [frm][3] - точки коллокации(центр ячейки)
    double** norm;// [frm][3] - нормаль к ячейке
    double*** root;// [frm][4][3] - ячейки
    double* s;//[frm] - площади ячеек


    int num_obj;// число объектов
    int num_mod;// число модулей 
    int num_frm;// число ячеек(общее)
    double max_diag;// шаг сетки

    int** beg_end;// [mod][2] номер 1 и посл ячейки на модуле с номером mod
    int** beg_end_obj;//[obj][2] 1 и посл модуль на каждом объекте с нoмером obj(сплошная сетка)

    void fill_TGrid(std::string filename); // функция заполнения информации о сетке

//------------------------------------------Constructor-------------------------------------------------------
    TGrid_DC_Full(std::string filename);//сквозная сетка на все объекты и все модули вместе

//----------------------------------------Copy Constructor----------------------------------------------------
    TGrid_DC_Full(const TGrid_DC_Full &obj);

//------------------------------------------Destructor--------------------------------------------------------
    ~TGrid_DC_Full();
};
#endif
