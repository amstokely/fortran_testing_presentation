
configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/config.hpp.in
        ${CMAKE_CURRENT_BINARY_DIR}/config.hpp
        @ONLY
)

add_library(mutual_information_cpp SHARED
        mutual_information.cpp
)

if (ENABLE_CUDA)
    target_sources(mutual_information_cpp PRIVATE
            mutual_information.cu
    )

    set_target_properties(mutual_information_cpp PROPERTIES
            CUDA_SEPARABLE_COMPILATION ON
    )
endif()

set_target_properties(mutual_information_cpp PROPERTIES
        LINKER_LANGUAGE CXX
)
target_include_directories(mutual_information_cpp PRIVATE
        ${CMAKE_CURRENT_BINARY_DIR}
)

target_include_directories(mutual_information_cpp INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<INSTALL_INTERFACE:include>
)


# ------------------------------------------------------------
# Fortran interface library
# ------------------------------------------------------------
if (ENABLE_CUDA)
    add_library(mutual_information SHARED
            mutual_information.f90
            cpp_mutual_information.f90
            cuda_mutual_information.f90
    )
else()
    add_library(mutual_information SHARED
            mutual_information.f90
            cpp_mutual_information.f90
            mutual_information_cuda_impl.cuh
    )
endif()

target_link_libraries(mutual_information PUBLIC
        mutual_information_cpp
)

target_include_directories(mutual_information INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<INSTALL_INTERFACE:include>
)
