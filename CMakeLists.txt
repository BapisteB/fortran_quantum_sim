cmake_minimum_required(VERSION 3.14)
project(my_proj LANGUAGES Fortran)

find_package(OpenMP REQUIRED COMPONENTS Fortran)
add_executable(my_prog main.f90 statevec.f90 gate.f90)

# Debug flags
target_compile_options(my_prog PRIVATE
    $<$<CONFIG:Debug>:-g -O0 -fcheck=all -fbacktrace -pg>
)

# Release flags
target_compile_options(my_prog PRIVATE
    $<$<CONFIG:Release>:-O3>
)

# Flags
set(FVERSION "-std=gnu")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FVERSION}")

# Standard lib
target_include_directories(my_prog PRIVATE /usr/local/include/fortran_stdlib/GNU-15.1.1)
target_link_libraries(my_prog PRIVATE /usr/local/lib/libfortran_stdlib.a)

# OpenMP
target_link_libraries(my_prog PRIVATE OpenMP::OpenMP_Fortran)
