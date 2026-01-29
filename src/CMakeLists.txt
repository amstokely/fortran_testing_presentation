add_library(boost_math STATIC boost_math.cpp)
set_target_properties(boost_math PROPERTIES
        LINKER_LANGUAGE CXX
)
add_library(mutual_information STATIC
        mutual_information.f90
        mutual_information_cpp.f90
)
target_include_directories(mutual_information INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<INSTALL_INTERFACE:include>
)