add_library(mutual_information_cuda SHARED mutual_information.cu)
set_target_properties(mutual_information_cuda PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        LINKER_LANGUAGE CUDA
)
target_include_directories(mutual_information_cuda INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
)

add_library(mutual_information_cpp SHARED mutual_information.cpp)
set_target_properties(mutual_information_cpp PROPERTIES
        LINKER_LANGUAGE CXX
)
target_include_directories(mutual_information_cpp INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
)
target_link_libraries(mutual_information_cpp PRIVATE mutual_information_cuda)

add_library(mutual_information SHARED
        mutual_information.f90
        cpp_mutual_information.f90
        cuda_mutual_information.f90
)
target_link_libraries(mutual_information PRIVATE mutual_information_cpp)
target_include_directories(mutual_information INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<INSTALL_INTERFACE:include>
)