cmake_minimum_required(VERSION 3.15)
project(numalz)

set(CMAKE_CXX_STANDARD 14)

find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
link_directories(${OpenCV_LIBRARY_DIRS} ${PROJECT_BINARY_DIR})

add_executable(numalz_matrix_test module_test/matrix_test.cpp Matrix.h)
add_executable(numalz_interop_test Interop.h module_test/interop_test.cpp )
add_executable(numalz_solver_test Solver.h module_test/solver_test.cpp)
add_executable(numalz_opencv_test module_test/opencv_test.cpp)
add_executable(numalz_axq module_test/q_2_2.cpp Matrix.h LinearAlgebra.h Solver.h Interop.h)

target_link_libraries(numalz_opencv_test ${OpenCV_LIBS})