cmake_minimum_required(VERSION 2.8)

#project(gatbl-tests)
include_directories( "../thirdparty/catch/single_include/catch2" )
add_executable(tests "main.cpp" "fastx.cpp" "sched.cpp")
#set_source_files_properties(sched.cpp PROPERTIES COMPILE_FLAGS -fsanitize=thread)
#set_target_properties(tests PROPERTIES LINK_FLAGS -fsanitize=thread )
target_link_libraries(tests  -lrt -lpthread)
