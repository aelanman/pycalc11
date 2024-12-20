cmake_minimum_required(VERSION 3.17.2...3.29)
project(${SKBUILD_PROJECT_NAME} LANGUAGES C Fortran)

#set(PYTHON_EXECUTABLE "/home/alanman/.venv/py311/bin/python")

find_package(
  Python 3.9
  COMPONENTS Interpreter Development.Embed Development.Module NumPy
  REQUIRED)

#set(CMAKE_CURRENT_BINARY_DIR "/home/alanman/temp/pycalc11/build")

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

set(f2py_module_name "calc11")
set(fortran_src_file "${CMAKE_BINARY_DIR}/calc11/src/calc11.f90")
set(generated_module_file ${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}${PYTHON_EXTENSION_MODULE_SUFFIX})

# Generate param11.in file in cmake build directory
file(READ ${CMAKE_SOURCE_DIR}/calc11/src/param11.i.in CONTENTS)
file(WRITE ${CMAKE_SOURCE_DIR}/calc11/src/param11.i ${CONTENTS})

add_subdirectory(calc11)

message(STATUS "SKBUILD_PROJECT_NAME=${SKBUILD_PROJECT_NAME}")
message(STATUS "F2PY_INCLUDE_DIR=${F2PY_INCLUDE_DIR}")
message(STATUS "CMAKE_CURRENT_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}")
message(STATUS "CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "generated_module_file=${generated_module_file}")
message(STATUS "fortran_source_file=${fortran_src_file}")
message(STATUS "f2py_module_name=${f2py_module_name}")
file(GLOB includes ${CMAKE_SOURCE_DIR}"/calc11/src/*.i")
#list(PREPEND includes ${generated_module_file})

#add_custom_target(${f2py_module_name} ALL
#  DEPENDS ${includes}
#  )

add_custom_command(
  OUTPUT calc11module.c calc11-f2pywrappers.f calc11-f2pywrappers2.f90
  DEPENDS ${includes}
  VERBATIM
  COMMAND "${Python_EXECUTABLE}" -m numpy.f2py
          "${fortran_src_file}" -m ${f2py_module_name} --lower
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
 
python_add_library(
        ${f2py_module_name}
       MODULE
       "${CMAKE_CURRENT_BINARY_DIR}/calc11-f2pywrappers2.f90"
       "${CMAKE_CURRENT_BINARY_DIR}/calc11-f2pywrappers.f"
       "${CMAKE_CURRENT_BINARY_DIR}/calc11module.c"
       "${fortran_src_file}"
       WITH_SOABI)
target_link_libraries(${f2py_module_name} PRIVATE fortranobject)

install(
    TARGETS ${f2py_module_name}
    DESTINATION ${SKBUILD_PROJECT_NAME}
)
