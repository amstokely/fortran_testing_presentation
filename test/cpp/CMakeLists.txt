add_executable(test_ksg_count_cpp test_ksg_count.cpp)
set_target_properties(
        test_ksg_count_cpp PROPERTIES
        LINKER_LANGUAGE CXX
)
target_link_libraries(
        test_ksg_count_cpp
        PRIVATE
        mutual_information_cpp
        test_utils_cpp
        project_options
)

add_test(
        NAME test_ksg_count_cpp
        COMMAND test_ksg_count_cpp
)