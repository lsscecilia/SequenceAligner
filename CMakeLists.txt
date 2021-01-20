cmake_minimum_required(VERSION 3.10)

#set project name
project(mapper VERSION "1.0")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(BINARY ${CMAKE_PROJECT_NAME})
set(SOURCE src/mapper.cpp src/alignment.cpp src/minimizer.cpp)
set(CONFIGURE_FILE src/config.h.in)
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})	
#set(PARSER_INCLUDES ${CMAKE_SOURCE_DIR}/3rdparty/bioparser)

#set executable
add_executable(${PROJECT_NAME} ${SOURCE}) 

configure_file( src/config.h.in ${CMAKE_BINARY_DIR}/bin/config.h)
include_directories( ${CMAKE_BINARY_DIR}/bin/ ) # Make sure it can be included...

if (NOT TARGET bioparser)
  add_subdirectory(3rdparty/bioparser EXCLUDE_FROM_ALL)
endif ()
#target_link_libraries(mapper PUBLIC bioparser)

add_library(
    alignment src/alignment.cpp
)

add_library(
    minimizer src/minimizer.cpp
)


target_link_libraries(
    	mapper
        bioparser
	alignment
	minimizer
)

#google test
option(mapper_tests "tests" ON)

if (mapper_tests)
	include_directories(${CMAKE_SOURCE_DIR}/src)
	add_executable(mapper_test test/alignmentTest.cpp)
	add_subdirectory(3rdparty/googletest EXCLUDE_FROM_ALL)
	target_link_libraries(mapper_test alignment minimizer gtest_main)
	
	
endif(mapper_tests)



