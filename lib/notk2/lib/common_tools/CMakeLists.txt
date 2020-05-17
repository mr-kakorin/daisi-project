cmake_minimum_required(VERSION 2.8)
project(common_tools)

set(COMM_TOOLS "COMM_TOOLS")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
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
add_definitions(-DBOOST_LOG_DYN_LINK=1)
find_package(Boost 1.55)

include_directories (include)
include_directories (lib/boostlog0/include)
#include_directories (include/common_tools)
include_directories (${Boost_INCLUDE_DIRS})

aux_source_directory(src SRC)

add_library(${PROJECT_NAME} SHARED ${SRC})
target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES} )

target_compile_definitions(${PROJECT_NAME} PRIVATE -DLOG_CHANNEL=COMM_TOOLS)
