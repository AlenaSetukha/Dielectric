cmake_minimum_required(VERSION 3.22)
project(Dielectric CXX)

set(CMAKE_CXX_STANDARD 11)

set(cmake_cxx_flags_debug "-Wall -Werror")
set(cmake_cxx_release "-Ofast -flto -march=native")

find_package(OpenMP REQUIRED)


set(DIEL_LIBRARY_SOURCE
    block.cpp
    common_type.cpp
    e0.cpp
    ED_Par.cpp
    element_geom.cpp
    epr_lib.cpp
    f_par.cpp
    get_area.cpp
    get_b.cpp
    get_einc.cpp
    get_epr.cpp
    get_field.cpp
    get_field_one_point.cpp
    get_hinc.cpp
    get_j.cpp
    Get_Matrix.cpp
    integral_par.cpp
    K.cpp
    kernel_lib.cpp
    Num_Par.cpp
    R.cpp
    )

add_library(field_library SHARED ${DIEL_LIBRARY_SOURCE})
target_link_libraries(field_library OpenMP::OpenMP_CXX)

add_executable(main main0_general.cpp)
target_link_libraries(main PUBLIC field_library)
