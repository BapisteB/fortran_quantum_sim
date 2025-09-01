cmake_minimum_required(VERSION 3.14)
project(my_proj LANGUAGES Fortran)

add_executable(my_prog main.f90 statevec.f90 gate.f90)

# Debug flags
target_compile_options(my_prog PRIVATE
    $<$<CONFIG:Debug>:-g -O0 -fcheck=all -fbacktrace>
)

# Release flags
target_compile_options(my_prog PRIVATE
    $<$<CONFIG:Release>:-O3>
)

target_include_directories(my_prog PRIVATE /usr/local/include/fortran_stdlib/GNU-15.1.1)
target_link_libraries(my_prog PRIVATE /usr/local/lib/libfortran_stdlib.a)
