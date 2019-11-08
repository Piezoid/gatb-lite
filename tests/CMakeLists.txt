list(APPEND test_sources "fastx.cpp")
list(APPEND test_sources "sched.cpp")
list(APPEND test_sources "kmer.cpp")
list(APPEND test_sources "memory.cpp")

add_executable(tests "main.cpp" ${test_sources})
set_source_files_properties(${test_sources} PROPERTIES COMPILE_FLAGS "-UNDEBUG -fsanitize=address -fsanitize=undefined -Werror -Wall -Wextra -Wconversion -pedantic-errors -Wno-gnu-zero-variadic-macro-arguments")
target_link_options(tests PRIVATE -fsanitize=address -fsanitize=undefined)
target_link_libraries(tests gatbl Catch2)

#set_source_files_properties(sched.cpp PROPERTIES COMPILE_FLAGS -fsanitize=thread)

file(GLOB_RECURSE test_fata_files "data/*")
add_custom_command(OUTPUT data
                   COMMAND echo "copying..."
                   COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/data data
                   DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/data ${test_fata_files})
add_custom_target(test_data DEPENDS data SOURCES ${test_fata_files})
add_dependencies(tests test_data)


include(CTest)
include(Catch)
catch_discover_tests(tests)
