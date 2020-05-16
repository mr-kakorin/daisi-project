cmake_minimum_required(VERSION 2.8)

aux_source_directory(accel/base SRC)
aux_source_directory(accel/Synchrotron SRC)
aux_source_directory(accel/RFQ SRC)
aux_source_directory(common SRC)
aux_source_directory(simulations SRC)

set(Boost_USE_STATIC_LIBS Off)
set(Boost_USE_MULTITHREAD On)
set(Boost_USE_STATIC_RUNTIME Off)

add_definitions(-DBOOST_LOG_DYN_LINK=1)

find_package(Boost 1.61 COMPONENTS serialization REQUIRED)
find_package(Armadillo REQUIRED)

include_directories (. ${Boost_INCLUDE_DIRS})
include_directories (. ${ARMADILLO_INCLUDE_DIRS})
include_directories (lib/include)
include_directories (common)

add_library(${PROJECT_NAME} SHARED ${SRC})

target_link_libraries(
    ${PROJECT_NAME} 
    ${Boost_LIBRARIES}
    ${ARMADILLO_LIBRARIES}
    ${VTK_LIBRARIES}
    notk)