find_package(PFUNIT REQUIRED)

add_pfunit_ctest(hashable_test
  TEST_SOURCES test_hashable.pf
  LINK_LIBRARIES history
)

add_pfunit_ctest(hash_table_test
  TEST_SOURCES test_hash_table.pf
  LINK_LIBRARIES history
)
