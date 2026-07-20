cmake_minimum_required(VERSION 3.16)

project(kmat VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(KMAT_BUILD_TESTS "Build kmat unit tests" ON)
option(KMAT_BUILD_BENCHES "Build kmat benchmark harness" ON)

include(FetchContent)
set(FETCHCONTENT_QUIET OFF)

FetchContent_Declare(
  CLI11
  GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
  GIT_TAG v2.4.2
  GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(CLI11)

set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
  Eigen3
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3.4.0
  GIT_SHALLOW TRUE
)
FetchContent_GetProperties(Eigen3)
if(NOT eigen3_POPULATED)
  FetchContent_Populate(Eigen3)
  if(NOT TARGET Eigen3::Eigen)
    add_library(Eigen3_Eigen INTERFACE)
    add_library(Eigen3::Eigen ALIAS Eigen3_Eigen)
    target_include_directories(Eigen3_Eigen SYSTEM INTERFACE "${eigen3_SOURCE_DIR}")
  endif()
endif()

if(KMAT_BUILD_TESTS)
  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.5.4
    GIT_SHALLOW TRUE
  )
  FetchContent_MakeAvailable(Catch2)
  include(CTest)
  include(Catch)
endif()

add_subdirectory(lib)
add_subdirectory(cli)

if(KMAT_BUILD_TESTS)
  add_subdirectory(tests)
endif()

if(KMAT_BUILD_BENCHES)
  if(NOT KMAT_BUILD_TESTS)
    include(CTest)
  endif()
  add_subdirectory(benches)
endif()
