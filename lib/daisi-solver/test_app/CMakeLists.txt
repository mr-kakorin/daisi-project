project(daisi_app)

set(Boost_USE_STATIC_LIBS Off)
set(Boost_USE_MULTITHREAD On)
set(Boost_USE_STATIC_RUNTIME Off)

add_definitions(-DBOOST_LOG_DYN_LINK=1)

add_definitions(-DDATA_PATH=${CMAKE_CURRENT_SOURCE_DIR}/data)

find_package(Boost 1.61 COMPONENTS log log_setup)

aux_source_directory( . SOURCES)

add_executable(${PROJECT_NAME} ${SOURCES})

target_link_libraries(${PROJECT_NAME}  ${Boost_LIBRARIES} daisi_solver pthread)