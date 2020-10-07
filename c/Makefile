

OBJS=svm_utils.o svm_fun.o stopwatch.o program_options.o svm_options.o svm_kernel_matrix.o \
        smo_solver.o svm_data.o svm_model.o svmtask.o

CFLAGS:=-O3 -ffast-math -mtune=native -march=native -funroll-loops -fopenmp

%.o: %.c
	$(CC) -g -c $(CFLAGS) $< -o $@

all:	$(OBJS)
	$(CC) $(CFLAGS)  -o svmtask $(OBJS) -lm

.phony: 
	all
	
clean:
	rm -fr *.o
	rm -fr svmtask	
