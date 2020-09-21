

OBJS=utils.o svm_options.o smo_solver.o svm_classify.o svmtask.o

CXXFLAGS:=-g -O3 -ffast-math -mavx2 -mtune=native -funroll-loops -std=c++17 -fopenmp

$(filter %.o,$(OBJS)): %.o: %.cpp
	$(CXX) -g -c $(CXXFLAGS) $< -o $@

all:	$(OBJS)
	rm -fr svmtask
	$(CXX) $(CXXFLAGS)  -o svmtask $(OBJS) -lm
	
clean:
	rm -fr *.o
	rm -fr svmtask	
