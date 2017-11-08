cmake_minimum_required(VERSION 2.8)

project(example1)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


# If the example is compiled independantly
if(NOT DEFINED GATBL_SOURCE_DIR)
    # An external project using gatb-lite would have this code

    set(GATBL_SOURCE_DIR "${CMAKE_SOURCE_DIR}/..")
    set(GATBL_INCLUDE_DIRECTORY "${GATBL_SOURCE_DIR}/include")
    include_directories(${GATBL_INCLUDE_DIRECTORY})

    find_package( ZLIB REQUIRED )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
endif()

add_executable(counter "counter.cpp")
target_link_libraries(counter -lpthread -lrt -lz)
