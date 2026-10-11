cmake_minimum_required(VERSION 3.16)
project(V8_9_Final_Scheme_C_Experiment LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Paper/benchmark build is intended for native single-process execution.
# Release mode uses strong CPU optimization; no threading library is linked.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(V89_COMMON_SOURCES
    src/TNR_2hop.cpp
    src/ET.cpp
)

add_executable(v89_final_benchmark
    main_v89_final_benchmark.cpp
    ${V89_COMMON_SOURCES}
)

target_compile_features(v89_final_benchmark PRIVATE cxx_std_17)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(v89_final_benchmark PRIVATE
        $<$<CONFIG:Release>:-O3>
        $<$<CONFIG:Release>:-march=native>
        $<$<CONFIG:Release>:-mtune=native>
        $<$<CONFIG:Release>:-DNDEBUG>
        $<$<CONFIG:Release>:-fomit-frame-pointer>
    )
    set_property(TARGET v89_final_benchmark PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
elseif(MSVC)
    target_compile_options(v89_final_benchmark PRIVATE
        $<$<CONFIG:Release>:/O2>
        $<$<CONFIG:Release>:/GL>
        $<$<CONFIG:Release>:/DNDEBUG>
    )
    set_property(TARGET v89_final_benchmark PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
endif()

option(BUILD_V89_CORRECTNESS_TEST "Build the validated V8.9 correctness test" ON)

if(BUILD_V89_CORRECTNESS_TEST)
    add_executable(test_v89_correctness
        test_v89_correctness.cpp
    )
    target_compile_features(test_v89_correctness PRIVATE cxx_std_17)
    target_include_directories(test_v89_correctness PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

    add_executable(test_tnr_end_to_end
        test_tnr_end_to_end.cpp
        ${V89_COMMON_SOURCES}
    )
    target_compile_features(test_tnr_end_to_end PRIVATE cxx_std_17)
    target_include_directories(test_tnr_end_to_end PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(test_v89_correctness PRIVATE
            $<$<CONFIG:Release>:-O3>
            $<$<CONFIG:Release>:-march=native>
            $<$<CONFIG:Release>:-mtune=native>
            $<$<CONFIG:Release>:-DNDEBUG>
        )
        target_compile_options(test_tnr_end_to_end PRIVATE
            $<$<CONFIG:Release>:-O2>
            $<$<CONFIG:Release>:-DNDEBUG>
        )
    endif()
endif()

message(STATUS "V8.9 final benchmark - Scheme C two-pointer Hub matching; clean query only")
message(STATUS "C++ standard: C++17")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Benchmark uses sequential single-threaded worker processes for per-case memory isolation")
