cmake_minimum_required(VERSION 3.6)
project(gui_client)

set(CMAKE_CXX_STANDARD 17)

if(CMAKE_BUILD_TYPE MATCHES Debug)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g3 -Wall -Wno-long-long -pedantic")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g0 -Os -march=native -flto")
endif()

add_subdirectory(../lib/daisi-solver "${CMAKE_CURRENT_BINARY_DIR}/daisi-solver")

add_definitions(-DBOOST_LOG_DYN_LINK)
add_definitions(-DBOOST_SP_USE_SPINLOCK -DBOOST_SP_USE_PTHREADS)
add_definitions(-DARMA_USE_CXX11 -DARMA_USE_LAPACK)
add_definitions(-DNUCL -DSIM)

include_directories(.)
include_directories(../lib/include)
include_directories(../lib)
set(BOOST_ROOT /usr/local)
set(Boost_USE_STATIC_LIBS Off)
set(Boost_USE_MULTITHREAD On)
set(Boost_USE_STATIC_RUNTIME Off)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_AUTOMOC ON)


find_package(Boost COMPONENTS serialization system REQUIRED)
find_package(Armadillo REQUIRED)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
find_package(OpenBLAS REQUIRED)

find_package(Qt5Widgets REQUIRED)
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

include_directories(${Boost_INCLUDE_DIRS})
include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${Qt5Widgets_INCLUDE_DIRS})
include_directories(${VTK_INCLUDE_DIRS})

include_directories(../lib/notk/lib/common_tools/include)
include_directories(../lib/daisi-solver/include)
include_directories(../lib/notk/include)
include_directories(common)

aux_source_directory(src SRC)

add_executable(${PROJECT_NAME}
        ${SRC}
        )

target_link_libraries(${PROJECT_NAME}
        ${Boost_LIBRARIES}
        ${ARMADILLO_LIBRARIES}
        ${VTK_LIBRARIES}
        Qt5::Widgets
        pthread
        ${OpenBLAS_LIB}
        daisi_solver
        )
