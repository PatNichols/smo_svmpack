

#include "utils.h"
#include "program_options.h"
#include "svm_options.h"
#include "smo_solver.h"
#include "svm_classify.h"

int main(int argc,char **argv)
{
    svm_options_t * options = svm_options_init(argc,argv);

    switch (options->task) {
    case 0:
        fprintf(stderr,"Training\n");
        smo_solver_train(options);
        break;
    case 1:
        fprintf(stderr,"Classify\n");
        svm_classify(options);
        break;
    case 2:
        fprintf(stderr,"Translating\n");
        svm_options_translate(options);
        break;
    default:
        fprintf(stderr,"No known task %d\n",options->task);
    }
    svm_options_free(options);
    return 0;
}