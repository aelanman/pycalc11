cmake_minimum_required(VERSION 3.15)

project(${SKBUILD_PROJECT_NAME}
  VERSION 1.0
  DESCRIPTION "CALC11 module"
  LANGUAGES C Fortran
  )


set(CMAKE_VERBOSE_MAKEFILE ON)

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)


#find_package(Python 3.8 REQUIRED)
#find_package(PythonExtensions REQUIRED)
#find_package(NumPy)
#find_package(F2PY REQUIRED)

execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c
          "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE)


# Generate param11.in file in cmake build directory
file(READ ${CMAKE_SOURCE_DIR}/calc11/src/param11.i.in CONTENTS)
file(WRITE ${CMAKE_SOURCE_DIR}/calc11/src/param11.i ${CONTENTS})


set(f2py_module_name "calc11")
set(fortran_src_file "${CMAKE_BINARY_DIR}/calc11/src/calc11.f90")
set(generated_module_file ${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}${PYTHON_EXTENSION_MODULE_SUFFIX})

add_subdirectory(calc11)

# get include files
file(GLOB includes ${CMAKE_SOURCE_DIR}"/calc11/src/*.i")
list(PREPEND includes ${generated_module_file})
add_custom_target(${f2py_module_name} ALL
  DEPENDS ${includes}
  )

add_library(fortranobject OBJECT "${F2PY_INCLUDE_DIR}/fortranobject.c")
target_link_libraries(fortranobject PUBLIC Python::NumPy)
target_include_directories(fortranobject PUBLIC "${F2PY_INCLUDE_DIR}")

add_custom_command(
  OUTPUT ${generated_module_file}
  COMMAND "${Python_EXECUTABLE}" -m numpy.f2py
    -m ${f2py_module_name}
    -c
    ${fortran_src_file}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

python_add_library(pycalc11 MODULE WITH_SOABI ${generated_module_file})

  #install(FILES ${generated_module_file} DESTINATION pycalc11)
