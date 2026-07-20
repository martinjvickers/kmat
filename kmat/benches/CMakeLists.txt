add_executable(kmat_bench bench_main.cpp)

target_link_libraries(kmat_bench
  PRIVATE
    kmat::kmat
)

target_compile_options(kmat_bench
  PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>
)

install(TARGETS kmat_bench RUNTIME DESTINATION bin)

add_test(
  NAME kmat_bench_laptop_tiny
  COMMAND kmat_bench --testdata "${CMAKE_SOURCE_DIR}/testdata" --panel tiny --profile laptop --threads 1
)
add_test(
  NAME kmat_bench_hpc_tiny
  COMMAND kmat_bench --testdata "${CMAKE_SOURCE_DIR}/testdata" --panel tiny --profile hpc --threads 2
)
set_tests_properties(kmat_bench_laptop_tiny kmat_bench_hpc_tiny PROPERTIES
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
)
