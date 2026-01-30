add_library(mutual_information_cpp STATIC mutual_information.cpp)
set_target_properties(mutual_information_cpp PROPERTIES
        LINKER_LANGUAGE CXX
)
target_include_directories(mutual_information_cpp INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
)
add_library(mutual_information STATIC
        mutual_information.f90
        cpp_mutual_information.f90
)
target_link_libraries(mutual_information PRIVATE mutual_information_cpp)
target_include_directories(mutual_information INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<INSTALL_INTERFACE:include>
)