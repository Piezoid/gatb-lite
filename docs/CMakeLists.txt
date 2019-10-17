find_program(SPHINX_EXECUTABLE NAMES sphinx-build
    HINTS
    $ENV{SPHINX_DIR}
    PATH_SUFFIXES bin
    DOC "Sphinx documentation generator"
)

if(Sphinx_FOUND)
    set(SPHINX_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}" CACHE PATH "Root of the output documentation")

    # configured documentation tools and intermediate build results
    set(SPHINX_BUILD_DIR "${SPHINX_OUTPUT_DIR}/_build" CACHE PATH "Intermediate output directory")

    # Sphinx cache with pickled ReST documents
    set(SPHINX_CACHE_DIR "${SPHINX_OUTPUT_DIR}/_doctrees" CACHE PATH "ReST cache directory")

    # HTML output directory
    set(SPHINX_HTML_DIR "${SPHINX_OUTPUT_DIR}/html" CACHE PATH "HTML output directory")

    mark_as_advanced(
      SPHINX_OUTPUT_DIR
      SPHINX_BUILD_DIR
      SPHINX_CACHE_DIR
      SPHINX_HTML_DIR
    )

    configure_file(
        "${SPHINX_INPUT_DIR}/conf.py.in"
        "${SPHINX_BUILD_DIR}/conf.py"
        @ONLY)

    add_custom_target("${PROJECT_NAME}_docs" ALL
        ${SPHINX_EXECUTABLE}
            -q -b html
            -c "${SPHINX_BUILD_DIR}"
            -d "${SPHINX_CACHE_DIR}"
            "${CMAKE_CURRENT_SOURCE_DIR}"
            "${SPHINX_HTML_DIR}"
        COMMENT "Building HTML documentation with Sphinx")
else()
    message(STATUS "Sphinx not found, the documentation will not be generated.")
endif()
