

#include "utils.h"
#include "program_options.h"
#include "svm_options.h"
#include "smo_solver.h"
#include "svm_classify.h"

int main(int argc,char **argv)
{
    int task;

    svm_options_t * options = svm_options_init(argc,argv);

    task = options->task;

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
        translate(options->data);
        break;
    default:
        fprintf(stderr,"No known task %d\n",task);
    }
    svm_options_free(options);
    return 0;
}