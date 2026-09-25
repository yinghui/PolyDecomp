# A project that uses an installed polydecomp:
#   cmake -S examples/consumer -B build-consumer -DCMAKE_PREFIX_PATH=<install prefix>
cmake_minimum_required(VERSION 3.20)
project(polydecomp_consumer LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
find_package(polydecomp 1.0 REQUIRED)
add_executable(consumer ../convex_and_tree.cpp)
target_link_libraries(consumer PRIVATE pd::polydecomp)
