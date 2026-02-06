cmake_minimum_required(VERSION 3.20)
project(fortran_testing_presentation LANGUAGES Fortran CXX C CUDA)
enable_testing()

find_package(CUDA REQUIRED)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")

set(PFUNIT_DIR "/home/astokely/software/pFUnit/build/installed/PFUNIT-4.10/cmake")
find_package(PFUNIT REQUIRED)

include_directories(external)
add_subdirectory(src)
add_subdirectory(test)

