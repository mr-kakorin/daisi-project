# Test framework
project(notk_test)

set(Boost_USE_STATIC_LIBS Off)
set(Boost_USE_MULTITHREAD On)
set(Boost_USE_STATIC_RUNTIME Off)

add_definitions(-DBOOST_LOG_DYN_LINK=1)

add_definitions(-DTEST_PATH=${CMAKE_CURRENT_SOURCE_DIR}/data)

find_package(Boost 1.55 COMPONENTS log log_setup)

find_package(GTest REQUIRED)
find_package(Threads REQUIRED)

include_directories(${GTEST_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(boostlogconf/include)

if(NOT TARGET boostlog_conf)
    add_subdirectory(boostlogconf)
endif()

aux_source_directory( . SOURCES)

add_executable(${PROJECT_NAME} ${SOURCES})
target_link_libraries(
                     ${PROJECT_NAME}
                     ${Boost_LIBRARIES}
                     ${GTEST_BOTH_LIBRARIES}  
                     ${CMAKE_THREAD_LIBS_INIT} 
                     notk
                     boostlog_conf
                     )

cotire(${PROJECT_NAME})
                   
GTEST_ADD_TESTS(${PROJECT_NAME} "" ${SOURCES})