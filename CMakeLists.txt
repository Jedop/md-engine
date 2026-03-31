cmake_minimum_required(VERSION 3.28)

project(MD_Engine)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB SRC_FILES src/*.cpp)

add_executable(MD_Engine ${SRC_FILES})
find_package(OpenMP REQUIRED)
target_link_libraries(MD_Engine PUBLIC OpenMP::OpenMP_CXX)
