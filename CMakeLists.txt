cmake_minimum_required(VERSION 2.8)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
set(CMAKE_USER_MAKE_RULES_OVERRIDE "${CMAKE_CURRENT_LIST_DIR}/cmake/overrides.cmake")

project(GATBL CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Sphinx)
find_package( ZLIB REQUIRED )

set(GATBL_INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/include)
include_directories( ${ZLIB_INCLUDE_DIRS} )
include_directories( ${GATBL_INCLUDE_DIRECTORY} )

add_subdirectory(tests)
add_subdirectory(benchs)
#add_subdirectory(examples)
