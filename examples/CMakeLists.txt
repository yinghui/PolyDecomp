# Examples; each is also run as a test. Standalone use against an installed
# polydecomp: see examples/consumer/CMakeLists.txt.
set(pd_examples double_with_retry convex_and_tree)
if(PD_WITH_GMP)
  list(APPEND pd_examples holes_to_convex slice_cheese)
endif()
foreach(ex ${pd_examples})
  add_executable(example_${ex} ${ex}.cpp)
  target_link_libraries(example_${ex} PRIVATE pd::polydecomp)
  pd_set_warnings(example_${ex})
  if(PD_BUILD_TESTS)
    add_test(NAME example_${ex} COMMAND example_${ex})
  endif()
endforeach()
