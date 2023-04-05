cmake_minimum_required(VERSION 3.24)
project(SIM2 C)

set(CMAKE_C_STANDARD 99)

add_executable(SIM2 main.c geom.c geom.h matriz.c matriz.h defs.h jacobi.c jacobi.h)
