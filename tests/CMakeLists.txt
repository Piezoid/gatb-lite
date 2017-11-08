cmake_minimum_required(VERSION 2.8)

#project(gatbl-tests)
include_directories( "../thirdparty/catch/single_include/catch2" )
add_executable(tests "main.cpp" "fastx.cpp")
target_link_libraries(tests  -lrt -lpthread)
