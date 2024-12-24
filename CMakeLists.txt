cmake_minimum_required(VERSION 3.17.2...3.29)
project(${SKBUILD_PROJECT_NAME} LANGUAGES C Fortran)

find_package(
  Python 3.8
  COMPONENTS Interpreter Development.Embed Development.Module NumPy
  REQUIRED)

set(CMAKE_VERBOSE_MAKEFILE ON)

# F2PY headers
execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c
          "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE)

add_library(fortranobject OBJECT "${F2PY_INCLUDE_DIR}/fortranobject.c")
target_link_libraries(fortranobject PUBLIC Python::NumPy)
target_include_directories(fortranobject PUBLIC "${F2PY_INCLUDE_DIR}")
set_property(TARGET fortranobject PROPERTY POSITION_INDEPENDENT_CODE ON)

set(fortran_src_file "${CMAKE_BINARY_DIR}/calc11/src/calc11.f90")

# Generate param11.in file in cmake build directory
file(READ ${CMAKE_SOURCE_DIR}/calc11/src/param11.i.in CONTENTS)
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/calc11/src/param11.i ${CONTENTS})

add_subdirectory(calc11)

file(GLOB includes ${CMAKE_SOURCE_DIR}"/calc11/src/*.i")

add_custom_command(
  OUTPUT calc11module.c calc11-f2pywrappers.f calc11-f2pywrappers2.f90
  DEPENDS ${includes}
  VERBATIM
  COMMAND "${Python_EXECUTABLE}" -m numpy.f2py
          "${fortran_src_file}" -m calc11 --lower
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

python_add_library(
       calc11
       MODULE
       "${CMAKE_CURRENT_BINARY_DIR}/calc11-f2pywrappers2.f90"
       "${CMAKE_CURRENT_BINARY_DIR}/calc11-f2pywrappers.f"
       "${CMAKE_CURRENT_BINARY_DIR}/calc11module.c"
       "${fortran_src_file}"
       WITH_SOABI)
target_link_libraries(calc11 PRIVATE fortranobject)

if ("${SKBUILD_STATE}" STREQUAL "editable")
        set(build_dest "${CMAKE_SOURCE_DIR}/pycalc11")
else()
        set(build_dest ${SKBUILD_PROJECT_NAME})
endif()

install(
    TARGETS calc11
    DESTINATION ${build_dest}
)
