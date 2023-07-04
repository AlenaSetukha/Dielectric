#ifndef _SOLVE_H_
#define _SOLVE_H_

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>
#include "element_geom.h"

/*=================================================

		  ROTATE MATRIX

================================================== */
template <typename P>
void decomp(int n, P** am, int* ip)
{
    int n1, k1, m;
    P t;
    n1 = n - 1;
    ip[n1] = 1;

    for (int k = 0; k < n1; k++)
    {
std::cout << k << std::endl;
        k1 = k + 1;
        m = k;
        for (int i = k1; i < n; i++)
        {
            if (abs_tmp(am[i][k]) > abs_tmp(am[m][k])) m = i;
        }

        ip[k] = m;
        if (m != k)   ip[n1] = -ip[n1];
        t = am[m][k];
        am[m][k] = am[k][k];
        am[k][k] = t;

        if (abs_tmp(t) == 0)
        {
            exit(1);
        }

        for (int i = k1; i < n; i++)
        {
            am[i][k] = -am[i][k] / t;
        }

        for (int j = k1; j < n; j++)
        {
            t = am[m][j];
            am[m][j] = am[k][j];
            am[k][j] = t;
            for (int i = k1; i < n; i++)
            {
                am[i][j] = am[i][j] + am[i][k] * t;
            }
        }
    }
}


/*=================================================

*                  MATRIX SOLVE

*=================================================*/

template <typename P>
void solve(int n, P** am, P* b, int* ip)
{
    int n1, k1, m, k2;
    P t;

    n1 = n - 1;

    for (int k = 0; k < n1; k++)
    {
        k1 = k + 1;
        m = ip[k];
        t = b[m];
        b[m] = b[k];
        b[k] = t;
        for (int i = k1; i < n; i++)
        {
	    b[i] = b[i] + am[i][k] * t;
        }
    }

    for (int kb = 1; kb <= n1; kb++)
    {
        k2 = n - kb;
        b[k2] = b[k2] / am[k2][k2];
        t = -b[k2];
        for (int i = 0; i < k2; i++)
        {
            b[i] = b[i] + am[i][k2] * t;
        }
    }
    b[0] = b[0] / am[0][0];
}

#endif
