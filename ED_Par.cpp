#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>
#include <string>

#include "ED_Par.h"
#include "constants.h"
#include "element_geom.h"

ED_Par::ED_Par(std::string filename)
{
    Constants c;
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::cout << "Read " << filename <<  " error" << std::endl;
        exit(1);
    }
    fin >> alpha >> omega >> alpha_start >> alpha_end;

    fin >> k_medium;

    m_d = new std::complex<double>[k_medium];
    eps_d = new std::complex<double>[k_medium];
    k = new std::complex<double>[k_medium];
    lambda = new std::complex<double>[k_medium];

    for (int i = 0; i < k_medium; i++)
    {
        fin >> eps_d[i];
        fin >> m_d[i];
        k[i] = omega * sqrt(eps_d[i] * m_d[i] * c.eps0 * c.m0);
        lambda[i] = c.pi * 2. / k[i];
    }

    k_vec[0] =  - abs_tmp(k[0]) * cos(alpha);
    k_vec[1] =  - abs_tmp(k[0]) * sin(alpha);
    k_vec[2] = 0.;
    fin.close();
}

