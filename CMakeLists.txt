cmake_minimum_required (VERSION 2.6)
SET(CMAKE_CXX_COMPILER "g++")
#find_package(openmp)
project (test3)
add_executable(test3 main.cpp autopara.cpp model3.cpp sdnsim.cpp attacker.cpp)
add_executable(autotest input_main.cpp autopara.cpp model3.cpp sdnsim.cpp attacker.cpp)

#target_link_libraries(test3 openmp)
#target_include_directories(test3 ./)
#add_library(test3 STATIC)
set(CMAKE_BUILD_TYPE Release)
include_directories(../eigen3/)
SET(CMAKE_CXX_FLAGS "-O0")
SET( CMAKE_CXX_FLAGS  " ${CMAKE_CXX_FLAGS} -std=gnu++0x")
