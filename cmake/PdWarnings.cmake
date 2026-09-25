# Warning flags shared by every polydecomp target.
function(pd_set_warnings target)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Wshadow)
  endif()
endfunction()
