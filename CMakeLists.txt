cmake_minimum_required(VERSION 2.8.11)
project(ZKrypt)

set(CMAKE_MACOSX_RPATH 1)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)


include_directories(include)

add_library(
	zkrypt
	SHARED

	src/bn254.c
	src/bn254_params.c
	src/plonk.c
)

SET_TARGET_PROPERTIES(zkrypt PROPERTIES VERSION 1.0 SOVERSION 1)

INSTALL(TARGETS zkrypt LIBRARY DESTINATION lib)


# tests
enable_testing()

add_executable(plonk_test tests/plonk_test.c)

target_link_libraries (plonk_test LINK_PUBLIC zkrypt)

target_link_libraries (plonk_test LINK_PUBLIC gmssl)

add_executable(bn254_test tests/plonk_test.c)

target_link_libraries (bn254_test LINK_PUBLIC zkrypt)

target_link_libraries (bn254_test LINK_PUBLIC gmssl)

add_test(NAME plonk_test		COMMAND plonk_test)

add_test(NAME bn254_test		COMMAND bn254_test)