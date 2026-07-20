add_executable(kmat_tests
  test_kmer.cpp
  test_matrix.cpp
  test_stats.cpp
  test_fastq.cpp
  test_count.cpp
  test_build_stream.cpp
  test_runtime.cpp
  test_integration.cpp
  test_integration_medium.cpp
)

target_link_libraries(kmat_tests
  PRIVATE
    kmat::kmat
    Catch2::Catch2WithMain
)

target_compile_definitions(kmat_tests
  PRIVATE
    KMAT_TESTDATA_DIR="${CMAKE_SOURCE_DIR}/testdata"
)

target_compile_options(kmat_tests
  PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>
)

catch_discover_tests(kmat_tests)

# Optional: run hidden [kmc] smoke when binaries are on PATH.
find_program(KMAT_KMC_EXECUTABLE kmc)
find_program(KMAT_KMC_TOOLS_EXECUTABLE kmc_tools)
if(KMAT_KMC_EXECUTABLE AND KMAT_KMC_TOOLS_EXECUTABLE)
  add_test(
    NAME kmat_kmc_count_smoke
    COMMAND kmat_tests "[kmc]"
  )
  set_tests_properties(kmat_kmc_count_smoke PROPERTIES
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    LABELS "kmc"
  )
endif()

# Shell integration smoke test (optional manual parity check).
add_test(
  NAME kmat_cli_pipeline
  COMMAND "${CMAKE_SOURCE_DIR}/tests/run_pipeline.sh"
          "${CMAKE_BINARY_DIR}/cli/kmat"
          "${CMAKE_SOURCE_DIR}/testdata"
)
set_tests_properties(kmat_cli_pipeline PROPERTIES WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")

add_test(
  NAME kmat_cli_pipeline_medium
  COMMAND "${CMAKE_SOURCE_DIR}/tests/run_pipeline_medium.sh"
          "${CMAKE_BINARY_DIR}/cli/kmat"
          "${CMAKE_SOURCE_DIR}/testdata/panel_k31_n72"
)
set_tests_properties(kmat_cli_pipeline_medium PROPERTIES WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
