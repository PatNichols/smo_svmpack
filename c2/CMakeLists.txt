cmake_minimum_required( VERSION 3.10 FATAL_ERROR )

project( svmpack VERSION 1.0)

# specify the C++ standard
set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED True)

add_executable( svmtask
	svm_utils.c
	svm_fun.c
	stopwatch.c
	program_options.c
	svm_options.c
	svm_kernel_matrix.c
        smo_solver.c
	svm_data.c
	svm_model.c
	svmtask.c)

target_compile_options( svmtask PUBLIC -O3 -ffast-math -march=native -fopenmp)
target_link_libraries( svmtask PUBLIC m)
target_link_options( svmtask PUBLIC -fopenmp)
	


