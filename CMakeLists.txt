cmake_minimum_required(VERSION 3.20)
project(vortex CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -O3")

add_subdirectory(libs)
add_subdirectory(fock)
add_subdirectory(test)
add_subdirectory(examples)
