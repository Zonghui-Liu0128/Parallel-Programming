cmake_minimum_required(VERSION 3.21)
project(WaveEquation)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

SET(CMAKE_CXX_COMPILER mpicxx)
SET(CMAKE_C_COMPILER  mpicc)

set(CMAKE_CXX_STANDARD 23)

add_executable(WaveEquation my_code.cpp)

