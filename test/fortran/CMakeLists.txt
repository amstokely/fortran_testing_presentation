add_pfunit_ctest(
        test_mutual_information
        TEST_SOURCES test_mutual_information.pf
        LINK_LIBRARIES test_utils mutual_information project_options
)

add_pfunit_ctest(
        test_ksg_counts
        TEST_SOURCES test_ksg_counts.pf
        LINK_LIBRARIES test_utils mutual_information project_options
)

add_pfunit_ctest(
        test_f90_ksg_count
        TEST_SOURCES test_f90_ksg_count.pf
        LINK_LIBRARIES test_utils mutual_information
)