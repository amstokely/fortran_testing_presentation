add_library(test_utils_cpp SHARED test_utils.cpp c_test_utils.cpp)
target_include_directories(test_utils_cpp INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<INSTALL_INTERFACE:include>
)
target_link_libraries(test_utils_cpp PUBLIC mutual_information_cpp)
set_target_properties(test_utils_cpp PROPERTIES
        LINKER_LANGUAGE CXX
)
add_library(test_utils SHARED test_utils.f90)
set_target_properties(test_utils PROPERTIES
        LINKER_LANGUAGE Fortran
)
target_include_directories(test_utils INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<INSTALL_INTERFACE:include>
)
target_link_libraries(test_utils PUBLIC test_utils_cpp mutual_information)
