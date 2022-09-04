cmake_minimum_required(VERSION 3.9)
project(isd C)
include_directories(${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/light_m4ri/include)

add_executable(isd
  light_m4ri/src/matrix.c
  src/bits.c
  src/dumer.c
  src/isd.c
  src/sort.c
  src/transpose.c
  src/xoroshiro128plus.c)

option(PGO "Use Profile-guided optimization (set this option to GEN, then run the executable, then recompile setting this option to USE)" OFF)

foreach(option
    "BENCHMARK"
    "DUMER_DOOM"
    "DUMER_DOOM"
    "DUMER_EPS"
    "DUMER_L"
    "DUMER_LUT"
    "DUMER_LUT_SHIFT"
    "DUMER_LW"
    "DUMER_LW"
    "DUMER_P")
  if(${option})
    target_compile_definitions(isd PUBLIC ${option}=${${option}})
  endif()
endforeach()

if (${CMAKE_C_COMPILER_ID} MATCHES "Clang" OR ${CMAKE_C_COMPILER_ID} STREQUAL "GNU")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Ofast -march=native -g3")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")
endif()

if (${CMAKE_C_COMPILER_ID} MATCHES "Clang")
  if(${PGO} STREQUAL "GEN")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-instr-generate=build/isd-%p.profraw")
  endif()
  if(${PGO} STREQUAL "USE")
    file(GLOB profraw_list "build/isd-*.profraw")
    execute_process(COMMAND llvm-profdata merge -output=build/isd.profdata ${profraw_list})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-instr-use=isd.profdata")
  endif()
elseif (${CMAKE_C_COMPILER_ID} STREQUAL "GNU")
  if(${PGO} STREQUAL "GEN")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-generate=build/isd.profdata")
  endif()
  if(${PGO} STREQUAL "USE")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-use=isd.profdata")
  endif()
endif()

set_target_properties(isd
  PROPERTIES
  C_STANDARD 11
  C_STANDARD_REQUIRED YES
  C_EXTENSIONS YES
  )

find_package(OpenMP)
if(OpenMP_C_FOUND)
  target_link_libraries(isd PUBLIC OpenMP::OpenMP_C)
endif()

include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_result)
if(ipo_result)
  set_target_properties(isd PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

find_library(MATH_LIBRARY m)
if(MATH_LIBRARY)
  target_link_libraries(isd PUBLIC ${MATH_LIBRARY})
endif()

file(GLOB_RECURSE ch_list *.c *.h)
add_custom_target(
  format
  COMMAND clang-format
  -i
  -verbose
  -style=google
  ${ch_list})
