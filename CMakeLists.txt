cmake_minimum_required (VERSION 2.6)
SET(CMAKE_CXX_COMPILER "g++")
#find_package(openmp)
project (test)
add_executable(test3 main.cpp model3.cpp sdnsim.cpp attacker.cpp)
#target_link_libraries(test3 openmp)
#target_include_directories(test3 ./)
#add_library(test3 STATIC)
set(CMAKE_BUILD_TYPE Release)
include_directories(../eigen3/)
SET(CMAKE_CXX_FLAGS "-O3")
SET( CMAKE_CXX_FLAGS  "-ftree-parallelize-loops=8 -fopenmp -fopenmp ${CMAKE_CXX_FLAGS} -std=gnu++0x")
        ## Optimize
  #  set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary

