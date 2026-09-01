cmake_minimum_required(VERSION 3.28)

project(MD_Engine LANGUAGES CXX CUDA)

set(CMAKE_CUDA_ARCHITECTURES 86)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB SRC_FILES src/*.cpp)

# CUDA source file
set(CUDA_SRC_FILES
    src/forces.cu
)

add_executable(MDEngine
    ${SRC_FILES}
    ${CUDA_SRC_FILES}
)

find_package(OpenMP REQUIRED)

target_link_libraries(MDEngine PUBLIC
    OpenMP::OpenMP_CXX
)

# Optimization for CUDA compilation
target_compile_options(MDEngine PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:-O3>
)


