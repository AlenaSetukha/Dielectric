#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>

#include "Num_Par.h"

Num_Par::Num_Par(const char* filename)
{
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::cout << "Read num_param.txt error" << std::endl;
        exit(1);
    }
    fin >> eps >> n_start >> n_start_seg >> p_max >> p_max_seg >> rs;
    fin.close();
}
