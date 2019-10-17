cmake_minimum_required(VERSION 3.5)

project(gatbl
    VERSION 0.0.1
    LANGUAGES CXX)
set(PROJECT_HOMEPAGE_URL "https://github.com/Piezoid/gatb-lite")
set(PROJECT_DESCRIPTION "Genomic Analysis ToolBox Lite, header only library for kmer based genomic analysis")
set(PROJECT_VENDOR "Maël Kerbiriou")
set(PROJECT_CONTACT "m431.kerbiriou@gmail.com")
string(COMPARE EQUAL ${CMAKE_SOURCE_DIR} ${PROJECT_SOURCE_DIR} SOURCE_IS_GATBL)
set(GATBL_STANDALONE ${SOURCE_IS_GATBL} CACHE BOOL "Build gatbl executables (tests and examples)")

include(GNUInstallDirs)

add_library(${PROJECT_NAME} INTERFACE) # Header only part, link to that to bring the includes in scope
target_include_directories(${PROJECT_NAME} INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(${PROJECT_NAME} INTERFACE pthread)
if(${CMAKE_VERSION} VERSION_GREATER "3.7.99")
    target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_11)
else()
    set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 11)
endif()
target_compile_options(${PROJECT_NAME} INTERFACE -msse4.2) #FIXME: multiple compilation unit + dipatch with vectorclass

install(DIRECTORY include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)



if(${GATBL_STANDALONE})
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
       add_compile_options (-fdiagnostics-color=always)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
       add_compile_options (-fcolor-diagnostics)
    endif ()

    include(CTest)
    option(BUILD_TESTING "Build tests and benchmarks" ON)
    add_subdirectory(thirdparty)

    if(${BUILD_TESTING})
        add_subdirectory(tests)
    endif()

    add_subdirectory(benchs)
    add_subdirectory(examples)
    add_subdirectory(docs)
endif()

