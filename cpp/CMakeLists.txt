cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

PROJECT( smopack VERSION 1.0 LANGUAGES CXX)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable( svmtask utils.cpp smo_solver.cpp svm_classify.cpp svm_options.cpp svmtask.cpp)

target_compile_options( svmtask PUBLIC -O3 -ffast-math -march=native -fopenmp)
target_link_libraries( svmtask PUBLIC m)
target_link_options( svmtask PUBLIC -fopenmp)
	







