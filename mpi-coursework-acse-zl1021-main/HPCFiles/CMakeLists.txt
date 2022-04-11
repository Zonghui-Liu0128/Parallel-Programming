cmake_minimum_required(VERSION 2.8)

# set the project name (and the name f)
project(WaveEquation)

# add the executable
add_executable(WaveEquation my_code.cpp)

find_package(MPI REQUIRED)

target_include_directories(WaveEquation PUBLIC ${MPI_INCLUDE_PATH})
target_link_libraries(WaveEquation ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(WaveEquation PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(WaveEquation PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()