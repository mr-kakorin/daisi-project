cmake_minimum_required(VERSION 3.12)
project(daisi_cli)

set(CMAKE_CXX_STANDARD 17)

if(CMAKE_BUILD_TYPE MATCHES Debug)
    if(MSVC)
        if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
            string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
        else()
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
        endif()
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g3 -Wall -Wno-long-long -pedantic")
    endif()
else()
    if(MSVC)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W0 /Ox")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g0 -march=native -Os -flto")
    endif()
endif()

#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fuse-ld=lld")

add_subdirectory(../lib/daisi-solver "${CMAKE_CURRENT_BINARY_DIR}/daisi-solver")

add_definitions(-DARMA_USE_CXX11 -DARMA_USE_LAPACK)
add_definitions(-DSIM)

include_directories(.)
include_directories(include)
include_directories(../lib/include)
include_directories(../lib)

#set(CMAKE_AUTOMOC ON)

add_definitions(-DBOOST_LOG_DYN_LINK)
add_definitions(-DBOOST_SP_USE_SPINLOCK -DBOOST_SP_USE_PTHREADS)
set(Boost_USE_STATIC_LIBS Off)
set(Boost_USE_MULTITHREAD On)
set(Boost_USE_STATIC_RUNTIME Off)

find_package(Boost 1.61 COMPONENTS serialization system log REQUIRED)
find_package(Armadillo REQUIRED)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE}) 

include_directories(${Boost_INCLUDE_DIRS})
include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${VTK_INCLUDE_DIRS})

include_directories(../lib/daisi-solver/include)
include_directories(src)

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