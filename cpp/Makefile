

SRCS=svmtask.cpp svm_classify.cpp smo_solver.cpp svm_options.cpp program_options.cpp utils.cpp
OBJS=svmtask.o svm_classify.o smo_solver.o svm_options.o program_options.o utils.o

CXXFLAGS:=-g -O3 -ffast-math -mavx2 -mtune=native -funroll-loops -std=c++17 -fopenmp

$(filter %.o,$(OBJS)): %.o: %.cpp
	$(CXX) -g -c $(CXXFLAGS) $< -o $@

all:	
	$(CXX) $(CXXFLAGS)  -c utils.cpp
	$(CXX) $(CXXFLAGS)  -c svm_options.cpp
	$(CXX) $(CXXFLAGS)  -c smo_solver.cpp
	$(CXX) $(CXXFLAGS)  -c svm_classify.cpp
	$(CXX) $(CXXFLAGS)  -c svmtask.cpp
	$(CXX) $(CXXFLAGS)  -o svmtask utils.o svm_options.o smo_solver.o svm_classify.o svmtask.o -lm
	
clean:
	rm -fr *.o
	rm -fr svmtask	