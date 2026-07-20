find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

add_library(kmat_lib
  src/kmer.cpp
  src/io.cpp
  src/log.cpp
  src/stats.cpp
  src/fasta.cpp
  src/fastq.cpp
  src/sequence.cpp
  src/presence.cpp
  src/count.cpp
  src/runtime.cpp
  src/table.cpp
  src/matrix.cpp
  src/pca.cpp
  src/gwas.cpp
  src/gene.cpp
  src/validate.cpp
)

add_library(kmat::kmat ALIAS kmat_lib)

target_include_directories(kmat_lib
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(kmat_lib
  PUBLIC
    Eigen3::Eigen
  PRIVATE
    ZLIB::ZLIB
    Threads::Threads
)

target_compile_options(kmat_lib
  PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>
)
