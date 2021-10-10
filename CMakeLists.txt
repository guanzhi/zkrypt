cmake_minimum_required(VERSION 2.8.11)
project(ZKrypt)

set(CMAKE_MACOSX_RPATH 1)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)


include_directories(include)

add_library(
	zkrypt
	SHARED

	lib/bn254.c
	lib/bn254_params.c
	lib/plonk.c
)

SET_TARGET_PROPERTIES(zkrypt PROPERTIES VERSION 1.0 SOVERSION 1)


