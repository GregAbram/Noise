cmake_minimum_required(VERSION 2.8.7)
project(Test)
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
include_directories(${Noise_SOURCE_DIR})

add_executable(sim sim.cpp)
add_executable(ramp ramp.cpp)
add_executable(vis vis.cpp)
add_executable(writer writer.cpp)


target_link_libraries(sim ${VTK_LIBRARIES} noise)
target_link_libraries(vis ${VTK_LIBRARIES})
target_link_libraries(ramp ${VTK_LIBRARIES})
target_link_libraries(writer ${VTK_LIBRARIES})

install(TARGETS sim vis ramp writer DESTINATION bin)

