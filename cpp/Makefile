

OBJS= utils.o svm_options.o smo_solver.o svm_classify.o svmtask.o

SRCS= utils.cpp svm_options.cpp smo_solver.cpp svm_classify.cpp svmtask.cpp

CXXFLAGS:=-g -O3 -ffast-math -mavx2 -mtune=native -funroll-loops -std=c++17 -fopenmp

%.o: %.cpp
	$(CXX) -g -c $(CXXFLAGS) $< -o $@

all:	$(OBJS)
	$(CXX) $(CXXFLAGS)  -o svmtask $(OBJS) -lm

.phony: all
	
clean:
	rm -fr *.o
	rm -fr svmtask	
