add_executable(kmat main.cpp)

target_link_libraries(kmat
  PRIVATE
    kmat::kmat
    CLI11::CLI11
)

target_compile_options(kmat
  PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>
)

install(TARGETS kmat RUNTIME DESTINATION bin)
