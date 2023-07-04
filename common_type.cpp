#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>

#include "element_geom.h"
#include "common_type.h"


void TGrid_DC_Full::fill_TGrid(std::string filename)
{
    max_diag = 0.;
    int t = 0, sum = 0;
    int kobj = 0, kmod = 0, sign = 0, nfm = 0, i11 = 0, i22 = 0;
    double diag1[3], diag2[3], len1, len2;

    std::ifstream fin2(filename);
    if (!fin2.is_open()) {
        std::cout << "Read geodat.dat error" << std::endl;
        exit(1);
    }
    fin2 >> kobj;

    for (int i = 0; i < kobj; i++) {
        fin2 >> kmod;
        beg_end_obj[i][0] = t;
        for (int j = 0; j < kmod; j++) {
            fin2 >> sign >> nfm >> i11 >> i22;
            beg_end[t][0] = sum;
            for (int k = 0; k < nfm; k++) {
                for (int g = 0; g < 4; g++) {
                    fin2 >> root[sum][g][0] >> root[sum][g][1] >> root[sum][g][2];
                }
                for (int g = 0; g < 3; g++)
                {
                    diag1[g] = root[sum][2][g] - root[sum][0][g];
                    diag2[g] = root[sum][3][g] - root[sum][1][g];
                }
                len1 = vec_length(diag1);
                len2 = vec_length(diag2);
                if (len1 > max_diag)
                {
                    max_diag = len1;
                }
                if (len2 > max_diag)
                {
                    max_diag = len2;
                }
                sum++;          //number of frames
            }
            beg_end[t][1] = sum - 1;
            t++;
        }
        beg_end_obj[i][1] = t - 1;
    }
    fin2.close();
//------------------------------------------s[], norm[][], rkt[][], tau[][][]---------------------------------
    double len, a[3], b[3], c[3], d[3], ac[3], bd[3];
    for (int i = 0; i < num_frm; i++) {
        for (int j = 0; j < 3; j++) {
            rkt[i][j] = (root[i][0][j] + root[i][1][j] +
                    root[i][2][j] + root[i][3][j]) * 0.25;
            for (int k = 0; k < 3; k++) {
                a[k] = (root[i][0][k] + root[i][1][k]) / 2.0;
                b[k] = (root[i][1][k] + root[i][2][k]) / 2.0;
                c[k] = (root[i][2][k] + root[i][3][k]) / 2.0;
                d[k] = (root[i][3][k] + root[i][0][k]) / 2.0;
                ac[k] = c[k] - a[k];
                bd[k] = d[k] - b[k];
            }

            vec_prod(ac, bd, norm[i]);
            len = vec_length(norm[i]);
            s[i] = len;
            for (int k = 0; k < 3; k++) {
                norm[i][k] /= len;
            }

            len = vec_length(ac);
            for (int k = 0; k < 3; k++) {
                tau[i][0][k] = ac[k] / len;
            }

            vec_prod(norm[i], tau[i][0], tau[i][1]);

            len = vec_length(tau[i][1]);
            for (int k = 0; k < 3; k++) {
                tau[i][1][k] /= len;
            }
        }
    }
}





//------------------------------------------Constructor-------------------------------------------------------
TGrid_DC_Full::TGrid_DC_Full(std::string filename)//сквозная сетка на все объекты и все модули вместе
{
    int smod = 0, sum = 0, kobj = 0, kmod = 0;
    int  sign, nfm, i11, i22;
    double x;

//------------------------Calculation of the number of frames--------------------
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::cout << "Read geodat.dat error" << std::endl;
        exit(-1);
    }

    fin >> kobj;
    for (int i = 0; i < kobj; i++) {
        fin >> kmod;
        smod += kmod;
        for (int j = 0; j < kmod; j++) {
            fin >> sign >> nfm >> i11 >> i22;//not used
            sum = sum + nfm;
            nfm = nfm * 12;
            for (int k = 0; k < nfm; k++) {
                fin >> x;
            }
        }
    }
    fin.close();
//-------------------------------Memory allocation-------------------------------
    num_obj = kobj;
    num_mod = smod;
    num_frm = sum;

    beg_end = new int*[num_mod];
    for (int i = 0; i < num_mod; i++) {
        beg_end[i] = new int[2];
    }

    beg_end_obj = new int*[num_obj];
    for (int i = 0; i < num_obj; i++) {
        beg_end_obj[i] = new int[2];
    }

    root = new double**[num_frm];
    for (int i = 0; i < num_frm; i++) {
        root[i] = new double*[4];
        for (int j = 0; j < 4; j++) {
            root[i][j] = new double[3];
        }
    }

    s = new double[num_frm];

    norm = new double*[num_frm];
    for (int i = 0; i < num_frm; i++) {
        norm[i] = new double[3];
    }

    rkt = new double*[num_frm];
    for (int i = 0; i < num_frm; i++) {
        rkt[i] = new double[3];
    }

    tau = new double**[num_frm];
    for (int i = 0; i < num_frm; i++) {
        tau[i] = new double*[2];
        for (int j = 0; j < 2; j++) {
            tau[i][j] = new double[3];
        }
    }
}


//----------------------------------------Copy Constructor----------------------------------------------------
TGrid_DC_Full::TGrid_DC_Full(const TGrid_DC_Full &obj)
{
    num_obj = obj.num_obj;
    num_mod = obj.num_mod;
    num_frm = obj.num_frm;
    max_diag = obj.max_diag;

    beg_end = new int*[num_mod];
    for (int i = 0; i < num_mod; i++) {
        beg_end[i] = new int[2];
        beg_end[i][0] = obj.beg_end[i][0];
        beg_end[i][1] = obj.beg_end[i][1];
    }

    beg_end_obj = new int*[num_obj];
    for (int i = 0; i < num_obj; i++) {
        beg_end_obj[i] = new int[2];
        beg_end_obj[i][0] = obj.beg_end_obj[i][0];
        beg_end_obj[i][1] = obj.beg_end_obj[i][1];
    }

    root = new double**[num_frm];
    for (int i = 0; i < num_frm; i++) {
        root[i] = new double*[4];
        for (int j = 0; j < 4; j++) {
            root[i][j] = new double[3];
            root[i][j][0] = obj.root[i][j][0];
            root[i][j][1] = obj.root[i][j][1];
            root[i][j][2] = obj.root[i][j][2];
        }
    }

    s = new double[num_frm];
    for (int i = 0; i < num_frm; i++)
    {
        s[i] = obj.s[i];
    }

    norm = new double*[num_frm];
    for (int i = 0; i < num_frm; i++) {
        norm[i] = new double[3];
        norm[i][0] = obj.norm[i][0];
        norm[i][1] = obj.norm[i][1];
        norm[i][2] = obj.norm[i][2];
    }

    rkt = new double*[num_frm];
    for (int i = 0; i < num_frm; i++) {
        rkt[i] = new double[3];
        rkt[i][0] = obj.rkt[i][0];
        rkt[i][1] = obj.rkt[i][1];
        rkt[i][2] = obj.rkt[i][2];
    }

    tau = new double**[num_frm];
    for (int i = 0; i < num_frm; i++) {
        tau[i] = new double*[2];
        for (int j = 0; j < 2; j++) {
            tau[i][j] = new double[3];
            tau[i][j][0] = obj.tau[i][j][0];
            tau[i][j][1] = obj.tau[i][j][1];
            tau[i][j][2] = obj.tau[i][j][2];
        }
    }
}


//------------------------------------------Destructor--------------------------------------------------------
TGrid_DC_Full::~TGrid_DC_Full()
{
    for (int i = 0; i < num_mod; i++) {
        delete[] beg_end[i];
    }
    delete[] beg_end;

    for (int i = 0; i < num_obj; i++) {
        delete[] beg_end_obj[i];
    }
    delete[] beg_end_obj;

    for (int i = 0; i < num_frm; i++) {
        for (int j = 0; j < 4; j++) {
            delete[] root[i][j];
        }
        delete[] root[i];
    }
    delete[] root;

    delete[] s;

    for (int i = 0; i < num_frm; i++) {
        delete[] norm[i];
    }
    delete[] norm;

    for (int i = 0; i < num_frm; i++) {
        delete[] rkt[i];
    }
    delete[] rkt;

    for (int i = 0; i < num_frm; i++) {
        for (int j = 0; j < 2; j++) {
            delete[] tau[i][j];
        }
        delete[] tau[i];
    }
    delete[] tau;
}
