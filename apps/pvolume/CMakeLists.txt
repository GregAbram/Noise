cmake_minimum_required(VERSION 2.8.7)
project(Test)
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
include_directories(${Noise_SOURCE_DIR})
add_executable(pvnoise pvolume.cxx)
target_link_libraries(pvnoise ${VTK_LIBRARIES} noise)
install(TARGETS pvnoise DESTINATION bin)

