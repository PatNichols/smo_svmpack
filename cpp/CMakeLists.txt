cmake_minimum_required(VERSION 3.0)

PROJECT( smopack VERSION 1.0 LANGUAGES CXX)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable( svmtask utils.cpp smo_solver.cpp svm_classify.cpp svm_options.cpp svmtask.cpp)





