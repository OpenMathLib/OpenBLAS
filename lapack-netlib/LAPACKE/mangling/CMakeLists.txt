cmake_minimum_required(VERSION 3.6)
project(MANGLING C Fortran)

add_executable(xintface Fintface.f Cintface.c)
