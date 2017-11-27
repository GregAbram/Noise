cmake_minimum_required(VERSION 2.8.7)
project(Test)
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
include_directories(${Noise_SOURCE_DIR})

add_executable(sim sim.cpp)
set_target_properties(sim PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
set_target_properties(sim PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
target_link_libraries(sim ${VTK_LIBRARIES} noise ${MPI_LIBRARIES})

add_executable(ramp ramp.cpp)
set_target_properties(ramp PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
set_target_properties(ramp PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
target_link_libraries(ramp ${VTK_LIBRARIES} ${MPI_LIBRARIES})

add_executable(vis vis.cpp)
set_target_properties(vis PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
set_target_properties(vis PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
target_link_libraries(vis ${VTK_LIBRARIES} ${MPI_LIBRARIES})

add_executable(writer writer.cpp)
set_target_properties(writer PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
set_target_properties(writer PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
target_link_libraries(writer ${VTK_LIBRARIES} ${MPI_LIBRARIES})

install(TARGETS sim vis ramp writer DESTINATION bin)

