cmake_minimum_required(VERSION 3.20)
project(polydecomp VERSION 1.4.3 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(PD_WITH_GMP    "Build the exact kernel (requires GMP/gmpxx)" ON)
option(PD_BUILD_TESTS "Build the test suite"                        ON)
option(PD_WITH_CGAL   "Build the CGAL baselines in pd_bench"        OFF)
option(PD_FUZZ        "Build libFuzzer targets (clang only)"        OFF)
option(PD_BUILD_EXAMPLES "Build the examples"                         ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(PdWarnings)

# ---------------------------------------------------------------- library
add_library(polydecomp
  src/version.cpp
  src/validate.cpp
  src/wkt.cpp
  src/obj.cpp
  src/leaves.cpp
  src/generate.cpp
)
add_library(pd::polydecomp ALIAS polydecomp)
target_include_directories(polydecomp PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_compile_definitions(polydecomp PUBLIC
  PD_VERSION_STRING="${PROJECT_VERSION}")
pd_set_warnings(polydecomp)

if(PD_WITH_GMP)
  include(PdGmp)
  target_link_libraries(polydecomp PUBLIC pd::gmp)
  target_compile_definitions(polydecomp PUBLIC PD_HAVE_GMP=1)
endif()

# ---------------------------------------------------------------- apps
foreach(app pd_validate pd_gen pd_decompose)
  add_executable(${app} apps/${app}.cpp)
  target_link_libraries(${app} PRIVATE pd::polydecomp)
  pd_set_warnings(${app})
endforeach()
if(PD_WITH_GMP)  # pd_slice needs the exact kernel
  add_executable(pd_slice apps/pd_slice.cpp)
  target_link_libraries(pd_slice PRIVATE pd::polydecomp)
  pd_set_warnings(pd_slice)
endif()

# pd_bench: our methods and the vendored PolyPartition, Clipper2 and Earcut baselines
# always, CGAL baselines with PD_WITH_CGAL (needs GMP). PolyPartition (MIT), Clipper2
# (Boost) and earcut.hpp (ISC) in third_party/ are used only by pd_bench and are not
# installed.
if(PD_WITH_GMP)
  add_library(pd_polypartition STATIC third_party/polypartition/polypartition.cpp)
  target_include_directories(pd_polypartition SYSTEM PUBLIC third_party/polypartition)
  add_library(pd_clipper2 STATIC
    third_party/clipper2/src/clipper.engine.cpp third_party/clipper2/src/clipper.offset.cpp
    third_party/clipper2/src/clipper.rectclip.cpp third_party/clipper2/src/clipper.triangulation.cpp)
  target_include_directories(pd_clipper2 SYSTEM PUBLIC third_party/clipper2/include)
  add_executable(pd_bench apps/pd_bench.cpp)
  target_link_libraries(pd_bench PRIVATE pd::polydecomp pd_polypartition pd_clipper2)
  target_include_directories(pd_bench SYSTEM PRIVATE third_party/earcut/include)
  target_compile_definitions(pd_bench PRIVATE PD_DATA_DIR="${PROJECT_SOURCE_DIR}/data")
  pd_set_warnings(pd_bench)
  if(PD_WITH_CGAL)
    find_package(CGAL REQUIRED)
    target_link_libraries(pd_bench PRIVATE CGAL::CGAL)
    target_compile_definitions(pd_bench PRIVATE PD_HAVE_CGAL=1)
  endif()
endif()

# ---------------------------------------------------------------- fuzzing
# Configure with the clang-fuzz preset (clang, sanitizers, PD_FUZZ=ON).
if(PD_FUZZ)
  if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(FATAL_ERROR "PD_FUZZ needs clang (libFuzzer)")
  endif()
  foreach(target fuzz_wkt fuzz_grid fuzz_slice)
    add_executable(${target} fuzz/${target}.cpp)
    target_link_libraries(${target} PRIVATE pd::polydecomp)
    target_compile_options(${target} PRIVATE -fsanitize=fuzzer)
    target_link_options(${target} PRIVATE -fsanitize=fuzzer)
  endforeach()
endif()

# ---------------------------------------------------------------- tests
if(PD_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
if(PD_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

# ---------------------------------------------------------------- docs
find_package(Doxygen QUIET)
if(DOXYGEN_FOUND)
  set(PD_HAVE_DOT NO)
  if(DOXYGEN_DOT_FOUND)
    set(PD_HAVE_DOT YES)
  endif()
  configure_file(docs/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(docs
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/docs
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMENT "API documentation in ${CMAKE_CURRENT_BINARY_DIR}/docs/html")
endif()

# ---------------------------------------------------------------- install
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
install(TARGETS polydecomp EXPORT polydecompTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY include/pd DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(TARGETS pd_validate pd_gen pd_decompose RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
if(PD_WITH_GMP)
  install(TARGETS pd_slice RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
install(EXPORT polydecompTargets NAMESPACE pd:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/polydecomp)
configure_package_config_file(cmake/polydecompConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/polydecompConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/polydecomp)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/polydecompConfigVersion.cmake
  COMPATIBILITY SameMajorVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/polydecompConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/polydecompConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/polydecomp)
install(FILES LICENSE README.md CHANGELOG.md CITATION.cff DESTINATION ${CMAKE_INSTALL_DOCDIR})
