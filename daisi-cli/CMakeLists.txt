cmake_minimum_required(VERSION 2.8)
project(daisi_cli)

set(CMAKE_CXX_STANDARD 17)

if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-std=c++17")
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
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g0 -march=native -Os")
    endif()
endif()

add_subdirectory(../lib/daisi-solver "../../solver-build-${CMAKE_BUILD_TYPE}")

add_definitions(-DBOOST_SP_USE_SPINLOCK -DBOOST_SP_USE_PTHREADS)
add_definitions(-DARMA_USE_CXX11 -DARMA_USE_LAPACK)
add_definitions(-DNUCL -DSIM)

include_directories(.)
include_directories(include)
include_directories(../lib/include)
include_directories(../lib)
set(BOOST_ROOT /usr/local)
set(Boost_USE_STATIC_LIBS Off)
set(Boost_USE_MULTITHREAD On)
set(Boost_USE_STATIC_RUNTIME Off)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_AUTOMOC ON) 

add_definitions(-DBOOST_LOG_DYN_LINK=1)


find_package(Boost 1.55 COMPONENTS serialization REQUIRED)
find_package(Armadillo REQUIRED)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE}) 

include_directories (${Boost_INCLUDE_DIRS})
include_directories (${ARMADILLO_INCLUDE_DIRS})
include_directories (${VTK_INCLUDE_DIRS})

include_directories (../lib/notk/lib/common_tools/include)
include_directories (../lib/daisi-solver/include)
include_directories (../lib/notk/include)
include_directories (src)

aux_source_directory(src SRC)

add_executable(${PROJECT_NAME} 
               ${SRC})

target_link_libraries(${PROJECT_NAME}
                      ${Boost_LIBRARIES}
                      ${ARMADILLO_LIBRARIES}
                      ${VTK_LIBRARIES}
                      pthread 
                      daisi_solver 
                     )