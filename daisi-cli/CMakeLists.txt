cmake_minimum_required(VERSION 3.12)
project(daisi_cli)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CMAKE_BUILD_TYPE MATCHES Debug)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g3 -Wall -Wno-long-long -pedantic")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g0 -flto -march=native -O3")
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fuse-ld=lld")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libc++ -lc++abi")
endif()

add_subdirectory(../lib/daisi-solver "${CMAKE_CURRENT_BINARY_DIR}/daisi-solver")

include_directories(.)
include_directories(include)
include_directories(../lib/include)
include_directories(../lib)

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
               ${SRC} src/ITasksCreator.cpp src/DiodeTasksCreator.cpp src/DiodeTasksCreator.h src/TriodeTasksCreator.cpp src/TriodeTasksCreator.h)
add_definitions(-D_GLIBCXX_USE_CXX11_ABI=1)
target_link_libraries(${PROJECT_NAME}
                      ${Boost_LIBRARIES}
                      ${ARMADILLO_LIBRARIES}
                      ${VTK_LIBRARIES}
                      pthread
                      daisi_solver
                     )