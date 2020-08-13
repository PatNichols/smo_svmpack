

SRCS=svmtask.c svm_classify.c smo_solver.c svm_options.c program_options.c utils.c
OBJS=svmtask.o svm_classify.o smo_solver.o svm_options.o program_options.o utils.o

CFLAGS:=-O3 -ffast-math -mavx2 -mtune=native -funroll-loops -fopenmp

$(filter %.o,$(OBJS)): %.o: %.c
	$(CC) -g -c $(CFLAGS) $< -o $@

all:	
	$(CC) $(CFLAGS)  -c utils.c
	$(CC) $(CFLAGS)  -c stopwatch.c
	$(CC) $(CFLAGS)  -c program_options.c	
	$(CC) $(CFLAGS)  -c svm_options.c
	$(CC) $(CFLAGS)  -c svm_kernel_matrix.c
	$(CC) $(CFLAGS)  -c smo_solver.c
	$(CC) $(CFLAGS)  -c svm_classify.c
	$(CC) $(CFLAGS)  -c svmtask.c
	$(CC) $(CFLAGS)  -o svmtask utils.o stopwatch.o program_options.o svm_options.o svm_kernel_matrix.o smo_solver.o svm_classify.o svmtask.o -lm
	
clean:
	rm -fr *.o
	rm -fr svmtask	