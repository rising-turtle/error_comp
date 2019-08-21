cmake_minimum_required(VERSION 2.8.3)
project(error_comp)

## System dependencies are found with CMake's conventions
# find_package(Boost REQUIRED COMPONENTS system)

set(CMAKE_CXX_FLAGS "-ggdb -g -O3 -w -std=c++11 -fpermissive") # -fPIC -fpermissive

find_package(Ceres REQUIRED)
# set(CERES_INCLUDE_DIRS "/home/davidz/work/3rdlib/ceres-solver/include")
include_directories(${CERES_INCLUDE_DIRS} "/usr/include/eigen3" "./src")

# link_libraries(${CERES_LIBRARIES})

add_subdirectory(src bin)