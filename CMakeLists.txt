cmake_minimum_required(VERSION 3.10)

#set project name
project(mapper VERSION 1.0)

set(CMAKE_CXX_STANDARD 11) 
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(SOURCE src/main.cpp)
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})	

#set executable
add_executable(${PROJECT_NAME} ${SOURCE})

if (NOT TARGET bioparser)
  add_subdirectory(<path_to_submodules>/bioparser EXCLUDE_FROM_ALL)
endif ()
target_link_libraries(<your_exe> bioparser)
