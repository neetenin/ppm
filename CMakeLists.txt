CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
project( ppm )
find_package( OpenCV REQUIRED )
add_executable( ppm ppm.cpp)
target_link_libraries( ppm ${OpenCV_LIBS} )
include_directories( ppm ${OpenCV_INCLUDE_DIRS})

set(CMAKE_CXX_COMPILER gcc-7)
set(CMAKE_CXX_FLAGS "-O3 -fopenmp -lstdc++ ${CMAKE_CXX_FLAGS}")