cmake_minimum_required(VERSION 3.12)
project(daisi_solver)

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libc++ -lc++abi")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fuse-ld=lld")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp=libomp")
elseif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()

set(CMAKE_CXX_STANDARD 17)

if(CMAKE_BUILD_TYPE MATCHES Debug)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g3 -Wall -Wno-long-long -pedantic")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g0 -Os -march=native -flto")
endif()

add_definitions(-DBOOST_LOG_DYN_LINK)
add_definitions(-DBOOST_SP_USE_SPINLOCK -DBOOST_SP_USE_PTHREADS)
set(Boost_USE_STATIC_LIBS Off)
set(Boost_USE_MULTITHREAD On)
set(Boost_USE_STATIC_RUNTIME Off)
add_definitions(-D_GLIBCXX_USE_CXX11_ABI=1)
add_definitions(-DARMA_USE_CXX11 -DARMA_USE_LAPACK)
add_definitions(-DSIM -DUSE_BLAS)

find_package(VTK REQUIRED)

include_directories (include)
include_directories (include/daisi-solver)
include_directories (lib/include)

include_directories (${VTK_INCLUDE_DIRS})

add_subdirectory(src)