#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>

#include "constants.h"
#include "common_type.h"
#include "element_geom.h"

int get_area(double* x, const TGrid_DC_Full& a)
{
    Constants c;
    double teta, vec[3];
    double sum = 0.;


    for (int i = 0; i < a.num_frm; i++)
    {
        vec[0] = a.root[i][0][0] - x[0];
        vec[1] = a.root[i][0][1] - x[1];
        vec[2] = a.root[i][0][2] - x[2];
//123
        teta = solid_angle(a.root[i][0], a.root[i][1], a.root[i][2], x);

        if (scal_prod(vec, a.norm[i]) > 0)
        {
            sum -= fabs(teta);
        } else {
            sum += fabs(teta);
        }

//134
        teta = solid_angle(a.root[i][0], a.root[i][2], a.root[i][3], x);

        if (scal_prod(vec, a.norm[i]) > 0)
        {
            sum -= fabs(teta);
        } else {
            sum += fabs(teta);
        }

    }

//Result
    if (sum > -0.5) {
        return 0;
    } else return -1;

}
