cmake_minimum_required(VERSION 3.20)
project(label_method_single)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE "Debug")
SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS}  -std=c++11 -Wall -g")
#SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")



add_executable(
        main
        src/main.cpp
        src/label_assign.cpp
        src/conversions.cpp
        src/ne_graph.cpp
)