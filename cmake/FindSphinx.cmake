find_program(SPHINX_EXECUTABLE NAMES sphinx-build
    HINTS
    $ENV{SPHINX_DIR}
    PATH_SUFFIXES bin
    DOC "Sphinx documentation generator"
)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Sphinx DEFAULT_MSG
    SPHINX_EXECUTABLE
)

if(NOT Sphinx_FOUND)
    message(STATUS "Sphinx not found, the documentation will not be generated.")
endif()

set(SPHINX_INPUT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/docs" CACHE PATH "Directory with the config.py.in file")

set(SPHINX_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/docs" CACHE PATH "Root of the output documentation")

# configured documentation tools and intermediate build results
set(SPHINX_BUILD_DIR "${SPHINX_OUTPUT_DIR}/_build" CACHE PATH "Intermediate output directory")

# Sphinx cache with pickled ReST documents
set(SPHINX_CACHE_DIR "${SPHINX_OUTPUT_DIR}/_doctrees" CACHE PATH "ReST cache directory")

# HTML output directory
set(SPHINX_HTML_DIR "${SPHINX_OUTPUT_DIR}/html" CACHE PATH "HTML output directory")

mark_as_advanced(
  SPHINX_INPUT_DIR
  SPHINX_OUTPUT_DIR
  SPHINX_BUILD_DIR
  SPHINX_CACHE_DIR
  SPHINX_HTML_DIR
  SPHINX_THEME
  SPHINX_THEME_DIR
)
