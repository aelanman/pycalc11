include(${CMAKE_CURRENT_SOURCE_DIR}/sourcelist.cmake)
include_directories(${calc11_SOURCE_DIR}/calc11/src)

#set(fortran_src_file "${CMAKE_CURRENT_BINARY_DIR}/calc11.f90")
#set(fortran_src_file "${CMAKE_BINARY_DIR}/calc11/src/calc11.f90")
#message(FATAL_ERROR ${fortran_src_file})
# Variable sources = list of .f source files
foreach (src_file ${sources})
  file(READ ${src_file} CONTENTS)
  file(APPEND ${fortran_src_file} "${CONTENTS}")
endforeach()

foreach (inc ${includes})
  configure_file(${inc} ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
endforeach()
