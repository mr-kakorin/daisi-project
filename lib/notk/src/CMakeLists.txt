cmake_minimum_required(VERSION 2.8)

project(notk)
aux_source_directory(. SRC)
aux_source_directory(./tools SRC)
aux_source_directory(./base SRC)
aux_source_directory(./solvers SRC)
aux_source_directory(./solvers/gaussearch SRC)
aux_source_directory(./solvers/uniform1d SRC)
aux_source_directory(./solvers/particleswarm SRC)
aux_source_directory(./solvers/boxwilson SRC)
aux_source_directory(./solvers/heavyball SRC)

set(Boost_USE_STATIC_LIBS Off)
set(Boost_USE_MULTITHREAD On)
set(Boost_USE_STATIC_RUNTIME Off)


find_package(Boost 1.55 COMPONENTS log_setup log system)

add_definitions(-DBOOST_LOG_DYN_LINK=1)

include_directories (
                    . 
                    ${Boost_INCLUDE_DIRS}
                    ${ARMADILLO_INCLUDE_DIRS}
                    )

add_library(${PROJECT_NAME} SHARED ${SRC})
target_link_libraries(
                     ${PROJECT_NAME}
                     ${Boost_LIBRARIES} 
                     common_tools
                     )
               
target_compile_definitions(${PROJECT_NAME} PRIVATE -DLOG_CHANNEL=NOTK)