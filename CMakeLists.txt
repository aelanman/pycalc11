cmake_minimum_required(VERSION 3.12)

project(calc11
  VERSION 1.0
  DESCRIPTION "CALC11 module"
  LANGUAGES C Fortran
  )

add_subdirectory(calc11)

find_package(Python 3.8 REQUIRED)
find_package(PythonExtensions REQUIRED)
find_package(NumPy)
find_package(F2PY REQUIRED)

set(f2py_module_name "calc11")
set(fortran_src_file "${CMAKE_SOURCE_DIR}/calc11/src/calc11.f90")
set(generated_module_file ${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}${PYTHON_EXTENSION_MODULE_SUFFIX})

add_custom_target(${f2py_module_name} ALL
  DEPENDS ${generated_module_file}
  )

add_custom_command(
  OUTPUT ${generated_module_file}
  COMMAND ${F2PY_EXECUTABLE}
    -m ${f2py_module_name}
    -c
    ${fortran_src_file}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )


install(FILES ${generated_module_file} DESTINATION pycalc11)
