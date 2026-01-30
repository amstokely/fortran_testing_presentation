cmake_minimum_required(VERSION 3.20)
project(fortran_testing_presentation LANGUAGES Fortran CXX C)
enable_testing()

set(CMAKE_CXX_STANDARD 23)

set(PFUNIT_DIR "/Users/astokely/software/pFUnit/build/installed/PFUNIT-4.10/cmake")
find_package(PFUNIT REQUIRED)

include_directories(external)
add_subdirectory(src)
add_subdirectory(test)

