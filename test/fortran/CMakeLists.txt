add_pfunit_ctest(
        test_mutual_information
        TEST_SOURCES test_mutual_information.pf
        LINK_LIBRARIES test_utils mutual_information
)
set_target_properties(test_mutual_information PROPERTIES
        LINKER_LANGUAGE Fortran
)


add_pfunit_ctest(
        test_ksg_count
        TEST_SOURCES test_ksg_count.pf
        LINK_LIBRARIES test_utils mutual_information
)