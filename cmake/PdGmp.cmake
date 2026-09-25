# Locate GMP and its C++ bindings (gmpxx) through pkg-config and expose them
# as the imported target pd::gmp.
find_package(PkgConfig REQUIRED)
pkg_check_modules(PD_GMPXX REQUIRED IMPORTED_TARGET gmpxx gmp)
if(NOT TARGET pd::gmp)
  add_library(pd::gmp ALIAS PkgConfig::PD_GMPXX)
endif()
message(STATUS "polydecomp: GMP ${PD_GMPXX_gmp_VERSION} / gmpxx ${PD_GMPXX_gmpxx_VERSION}")
