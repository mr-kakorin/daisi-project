# Test framework

project(notk_test)

find_package(Boost 1.55 COMPONENTS thread system REQUIRED)

add_definitions(-DBOOST_LOG_DYN_LINK=1)

add_definitions(-DTEST_PATH=${CMAKE_CURRENT_SOURCE_DIR}/data)

find_package(GTest REQUIRED)
find_package(Threads REQUIRED)

include_directories(${GTEST_INCLUDE_DIRS})
include_directories(boostlogconf/include)
add_definitions(-DBOOST_LOG_DYN_LINK=1)
aux_source_directory( . SOURCES)

add_executable(${PROJECT_NAME} ${SOURCES})
target_link_libraries(
                     ${PROJECT_NAME}
                     ${Boost_LIBRARIES}
                     ${GTEST_BOTH_LIBRARIES}  
                     ${CMAKE_THREAD_LIBS_INIT} 
                     notk
                     boostlog0
                     )

GTEST_DISCOVER_TESTS(${PROJECT_NAME} "" ${SOURCES})
