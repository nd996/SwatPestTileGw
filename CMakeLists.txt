cmake_minimum_required(VERSION 3.1)

project(SWAT_681 VERSION 0.1
        DESCRIPTION "SWAT Fortran program"
        LANGUAGES Fortran)
enable_language(Fortran)

# Currently setting the Fortran compiler to use -std=gnu, change this if you
# want a specific standard
set(FVERSION "-std=gnu")
set(CMAKE_Fortran_FLAGS "-ffree-line-length-0 -ffixed-line-length-0 -g -O2 ${CMAKE_Fortran_FLAGS} ${FVERSION}")
set(CMAKE_Fortran_COMPILER "gfortran")

# Source code
add_subdirectory(src)
