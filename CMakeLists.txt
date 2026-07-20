cmake_minimum_required(VERSION 3.15)
project(CH_TNR LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "No build type selected, defaulting to Release")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
endif()

option(USE_PGO_GENERATE "Build with PGO instrumentation to generate profile data" OFF)
option(USE_PGO_USE "Build using existing PGO profile data" OFF)

set(BASE_CXX_FLAGS "-O3 -DNDEBUG -march=native -mtune=native -fno-stack-protector -fno-exceptions -fno-rtti -funroll-loops")

# 注入 PGO 标志
if(USE_PGO_GENERATE)
    message(STATUS "🚀 [PGO Mode] GENERATE: Instrumentation enabled.")
    set(CMAKE_CXX_FLAGS_RELEASE "${BASE_CXX_FLAGS} -fprofile-generate")
elseif(USE_PGO_USE)
    message(STATUS "🚀 [PGO Mode] USE: Applying profile data for extreme optimization.")
    set(CMAKE_CXX_FLAGS_RELEASE "${BASE_CXX_FLAGS} -fprofile-use -fprofile-correction")
else()
    message(STATUS "✅ [Standard Mode] Release build with static optimizations.")
    set(CMAKE_CXX_FLAGS_RELEASE "${BASE_CXX_FLAGS}")
endif()
set(BASE_LINKER_FLAGS "-fuse-ld=lld -Wl,--icf=all")

if(USE_PGO_GENERATE)
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${BASE_LINKER_FLAGS} -fprofile-generate")
elseif(USE_PGO_USE)
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${BASE_LINKER_FLAGS} -fprofile-use")
else()
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${BASE_LINKER_FLAGS}")
endif()

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)

add_executable(
        CH_TNR
        TNR_2hop.cpp
        ET.cpp
        TNR_2hop.h
        Fast2Hop.h
        ET.h
        Graph.h
        main.cpp

)
