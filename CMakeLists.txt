cmake_minimum_required(VERSION 3.12)

project(calc11
  VERSION 1.0
  DESCRIPTION "CALC11 module"
  LANGUAGES C Fortran
  )


find_package(Python 3.8 REQUIRED)
find_package(PythonExtensions REQUIRED)
find_package(NumPy)
find_package(F2PY REQUIRED)

#set(DATADIR "")
#set(DE421_FILE "")

# Generate param11.in file in cmake build directory
file(READ ${CMAKE_SOURCE_DIR}/calc11/src/param11.i.in CONTENTS)
#string(REPLACE "@prefix@/share/difxcalc" "${DATADIR}" CONTENTS ${CONTENTS})
#string(REPLACE "@de421_file@" "${DATADIR}" CONTENTS ${CONTENTS})
file(WRITE ${CMAKE_SOURCE_DIR}/calc11/src/param11.i ${CONTENTS})
#file(WRITE ${CMAKE_SOURCE_DIR}/test.in ${CONTENTS})


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

add_custom_command(
  OUTPUT ${generated_module_file}
  COMMAND ${F2PY_EXECUTABLE}
    -m ${f2py_module_name}
    -c
    ${fortran_src_file}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )


install(FILES ${generated_module_file} DESTINATION pycalc11)
