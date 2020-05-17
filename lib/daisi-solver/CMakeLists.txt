cmake_minimum_required(VERSION 3.12)
project(daisi_solver)

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp=libomp")
    find_package(OpenMP)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
elseif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()

set(CMAKE_CXX_STANDARD 17)

#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fuse-ld=lld")

if(CMAKE_BUILD_TYPE MATCHES Debug)
    if(MSVC)
        if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
            string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
        else()
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
        endif()
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g3 -Wall -Wno-long-long -pedantic")
    endif()
else()
    if(MSVC)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W0 /Ox")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g0 -Os -march=native -flto")
    endif()
endif()

add_definitions(-DBOOST_LOG_DYN_LINK)
add_definitions(-DBOOST_SP_USE_SPINLOCK -DBOOST_SP_USE_PTHREADS)
set(Boost_USE_STATIC_LIBS Off)
set(Boost_USE_MULTITHREAD On)
set(Boost_USE_STATIC_RUNTIME Off)

add_definitions(-DARMA_USE_CXX11 -DARMA_USE_LAPACK)
add_definitions(-DSIM -DUSE_BLAS)

find_package(VTK REQUIRED)

include_directories (include)
include_directories (include/daisi-solver)
include_directories (lib/include)

include_directories (${VTK_INCLUDE_DIRS})

add_subdirectory(src)