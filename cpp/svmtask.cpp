#include <omp.h>
#include "utils.hpp"
#include "program_options.hpp"
#include "svm_options.hpp"
#include "smo_solver.hpp"
#include "svm_classify.hpp"



int main(int argc,char **argv)
{
    int task;

    svm_options options(argc,argv);
 
    task = options.task;

    switch (task) {
    case 0:
        fprintf(stderr,"Training\n");
        smo_solver_train(options);
        break;
    case 1:
        fprintf(stderr,"Classify\n");
        svm_classify(options);
        break;
    case 2:
        options.translate();
        break;
    default:
        fprintf(stderr,"No known task %d\n",task);
    }
    return 0;
}