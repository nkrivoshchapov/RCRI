cmake_minimum_required(VERSION 3.0)
enable_language(Fortran)
project(libtlc)

set(CMAKE_CXX_STANDARD 14)
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
find_library(GFORTRAN_LIBRARY NAMES gfortran)
add_library(tlc SHARED tlc_src/main.cpp tlc_src/sturm.c tlc_src/tripep_closure.f90 tlc_src/main_pep.f90)
target_link_libraries(tlc gfortran)

install(TARGETS tlc
        LIBRARY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/rcrilib/Solvers/)