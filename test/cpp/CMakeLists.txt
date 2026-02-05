add_executable(test_ksg_count_cpp test_ksg_count.cpp)
set_target_properties(
    test_ksg_count_cpp PROPERTIES
    LINKER_LANGUAGE CXX
)
target_link_libraries(test_ksg_count_cpp PRIVATE test_utils_cpp mutual_information_cpp mutual_information_cuda)
add_test(
    NAME test_ksg_count_cpp
    COMMAND test_ksg_count_cpp
)