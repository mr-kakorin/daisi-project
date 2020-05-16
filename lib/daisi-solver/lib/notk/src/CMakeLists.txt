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

add_library(${PROJECT_NAME} SHARED ${SRC})
target_link_libraries(
                     ${PROJECT_NAME}
                     ${Boost_LIBRARIES} 
                     common_tools
                     )
               
target_compile_definitions(${PROJECT_NAME} PRIVATE -DLOG_CHANNEL=NOTK)
