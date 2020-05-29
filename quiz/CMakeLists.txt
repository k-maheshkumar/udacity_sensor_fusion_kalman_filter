cmake_minimum_required (VERSION 3.0)
project (kalman_example)
 
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
 
add_executable (kalman_1d src/kalman_1d/kalman.cpp)
target_link_libraries (kalman_1d Eigen3::Eigen)


