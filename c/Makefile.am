bin_PROGRAMS = svmtask

svmtask_SOURCES = svm_utils.c svm_fun.c stopwatch.c program_options.c svm_options.c svm_kernel_matrix.c \
        smo_solver.c svm_data.c svm_model.c svmtask.c

AM_CFLAGS =
svmtask_CFLAGS = -O3 -ffast-math -march=native -fopenmp
svmtask_LFLAGS = -fopenmp 
svmtask_LDADD  = -lm
