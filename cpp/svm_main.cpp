#include <cstdlib>
#include <iostream>
#include <string>
#include "svm_data.hpp"
#include "svm_options.hpp"
#include "smo_solver.hpp"
#include "svm_classify.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace svmpack {

void adjust_threads(int nth)
{
#ifdef _OPENMP
    int nth_ = omp_num_threads();
    if (nth_ != nth) omp_set_num_threads(nth);
#else
    if ( nth!=0) {
        std::cout << "open mp is not supported\n";
    }
#endif
}

}

int main(int argc,char **argv)
{
    try {
        svmpack::svm_options options(argc,argv);
        int task = options.task;
//        adjust_threads(options.nths);
        switch (task)
        {
        case 0:
            // train
            std::cout << "training\n";
            svmpack::smo_solver_train(options);
            break;
        case 1:
            // validate
            std::cout << "validating\n";
            svmpack::svm_validate(options);
            break;
        case 2:
            std::cout << "classifying\n";
            svmpack::svm_classify(options);
            break;
        case 3:
            std::cout << "translating\n";
            svmpack::svm_translate(options.data);
            break;
        default:
            std::cerr << "unknown task " << task << "\n";
            exit(-1);
        }
    }
    catch (std::exception& e)
    {
        std::cerr << e.what() << "\n";
        exit(-1);
    }
    return EXIT_SUCCESS;
}