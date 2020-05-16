cmake_minimum_required(VERSION 2.8)
project(notk)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    add_compile_options(-std=c++14)
endif()

if(CMAKE_BUILD_TYPE MATCHES Debug)
    if(MSVC)
        if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
            string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
        else()
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
        endif()
    elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g3 -Wall -Wno-long-long -pedantic")
    endif()
else()
    if(MSVC)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W0 /Ox")
    elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -flto")
    endif()
endif()
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake")

add_definitions(-DBOOST_LOG_DYN_LINK)
find_package(Boost 1.55 COMPONENTS log_setup log system)

add_definitions(-DBOOST_MPL_CFG_NO_PREPROCESSED_HEADERS)
add_definitions(-DBOOST_MPL_LIMIT_VECTOR_SIZE=30)

if(NOT TARGET common_tools)
add_subdirectory(lib/common_tools)
endif()

include_directories (include)
include_directories (lib/boostlog0/include)
include_directories (lib/common_tools/include)
include_directories (lib/serreflection/include)
include_directories (${Boost_INCLUDE_DIRS})


add_subdirectory(src)

#option(BUILD_TEST "Build and run test" ON)

#if(BUILD_TEST)
#    enable_testing()
#    
#    if(NOT TARGET boostlog0)
#        add_subdirectory(lib/boostlog0)
#    endif()
#
#    add_subdirectory(test)
#endif()
