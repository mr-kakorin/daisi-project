cmake_minimum_required(VERSION 2.8)
project(serreflection.tests)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++11" COMPILER_SUPPORTS_CXX11)

if (NOT COMPILER_SUPPORTS_CXX11)
    message(FATAL_ERROR "${PROJECT_NAME} build required c++11 support (gcc 4.7+). Please use different compiler")
endif()

if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fPIC")
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
option(BUILD_TESTS "Link and run test" ON)

set(Boost_USE_STATIC_LIBS Off)
set(Boost_USE_MULTITHREAD On)
set(Boost_USE_STATIC_RUNTIME Off)

add_definitions(-DBOOST_LOG_DYN_LINK=1)

find_package(Boost 1.55 COMPONENTS log log_setup REQUIRED)

include_directories (. ${Boost_INCLUDE_DIRS})

find_package(GTest REQUIRED)
find_package(Threads REQUIRED)

aux_source_directory(. SRCTEST)

if(NOT TARGET boostlog_conf)
    add_subdirectory(boostlogconf)
endif()

include_directories (boostlogconf/include)

add_executable(${PROJECT_NAME} ${SRCTEST})
target_link_libraries(
                      ${PROJECT_NAME}
                      ${CMAKE_THREAD_LIBS_INIT}
                      ${Boost_LIBRARIES}
                      ${GTEST_BOTH_LIBRARIES}
                      boostlog_conf
                      )
GTEST_ADD_TESTS(${PROJECT_NAME} "" ${SRCTEST})

target_compile_definitions(${PROJECT_NAME} PRIVATE -DLOG_CHANNEL=TEST_SRFL)