# Add StateSetMPITest to the test suite
ecbuild_add_test( TARGET test_base_stateset_mpi
                  SOURCES StateSetMPITest.cc
                  ARGS "--gtest_output=xml:test_base_stateset_mpi.xml"
                  LIBS oops
                  TEST_DEPENDS test_base
                  MPI 4 )