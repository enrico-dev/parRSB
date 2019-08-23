cmake_minimum_required(VERSION 3.4 FATAL_ERROR)
project(parRSB C)

option(GENMAP_PAUL "Use Paul Fischer's genmap for parRSB" ON)

add_library(parRSB
    src/genmap.c
    src/genmap-vector.c
    src/genmap-handle.c
    src/genmap-comm.c
    src/genmap-eigen.c
    src/genmap-laplacian.c
    src/genmap-lanczos.c
    src/genmap-rsb.c
    src/parrsb-binsort.c
    src/parrsb-histogram.c
    src/genmap-chelpers.c
    src/parRSB.c )

add_subdirectory(3rd_party/gslib/gslib)

target_include_directories(parRSB PRIVATE
    src
    3rd_party/gslib/gslib/src
    ${CMAKE_CURRENT_BINARY_DIR}/3rd_party/gslib/gslib)
target_link_libraries(parRSB PRIVATE gs)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  target_compile_definitions(parRSB PRIVATE -DGENMAP_DEBUG)
endif()

if(${GENMAP_PAUL})
  target_compile_definitions(parRSB PRIVATE -DGENMAP_PAUL)
endif()
