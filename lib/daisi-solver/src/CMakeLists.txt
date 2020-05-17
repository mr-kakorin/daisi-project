cmake_minimum_required(VERSION 3.12)

aux_source_directory(common SRC)
aux_source_directory(simulations SRC)

find_package(Boost 1.61 COMPONENTS serialization system log REQUIRED)
find_package(Armadillo REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})
include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(../lib/include)
include_directories(../include)
include_directories(common)

add_library(${PROJECT_NAME} SHARED ${SRC})

target_link_libraries(
    ${PROJECT_NAME} 
    ${Boost_LIBRARIES}
    ${ARMADILLO_LIBRARIES}
    ${VTK_LIBRARIES}
    )