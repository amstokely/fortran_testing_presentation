cmake_minimum_required(VERSION 3.20)
project(fortran_testing_presentation LANGUAGES Fortran CXX C)
set(CMAKE_CXX_STANDARD 23)
enable_testing()

set(PFUNIT_DIR CACHE PATH "Path to PFUNIT installation")
find_package(PFUNIT REQUIRED)

option(ENABLE_CUDA "Enable CUDA support" OFF)

add_library(project_options INTERFACE)

target_compile_options(project_options INTERFACE
        $<$<COMPILE_LANGUAGE:Fortran>:-cpp>
        $<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-none>
)
if (ENABLE_CUDA)
    enable_language(CUDA)
    find_package(CUDA REQUIRED)
    target_compile_definitions(project_options INTERFACE CUDA_SUPPORT)
endif ()


include_directories(external)

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(slides)
