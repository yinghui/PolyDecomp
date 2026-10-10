add_library(pd_doctest_main STATIC test_main.cpp)
target_include_directories(pd_doctest_main PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/third_party)

add_executable(pd_tests
  test_smoke.cpp
  test_kernel.cpp
  test_geometry.cpp
  test_validate.cpp
  test_io.cpp
  test_verify.cpp
  test_stage2.cpp
  test_stage1.cpp
  test_generate.cpp
  test_pipeline.cpp
  test_differential.cpp
  test_symbolic.cpp
)
target_link_libraries(pd_tests PRIVATE pd::polydecomp pd_doctest_main)
target_compile_definitions(pd_tests PRIVATE PD_DATA_DIR="${PROJECT_SOURCE_DIR}/data")
pd_set_warnings(pd_tests)

add_test(NAME pd_tests COMMAND pd_tests)
