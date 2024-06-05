#ifndef SVM_OPTIONS_H
#define SVM_OPTIONS_H
#include <cstdlib>
#include <string>
#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "putils_program_options.hpp"

namespace svmpack {
struct svm_options {
    double eps;
    double cost;
    double kc1;
    double kc2;
    int kpow;
    int ktype;
    int kscale;
    int cache_size;
    int cache_mem;
    int cache_precompute;
    int max_its;
    int task;
    int nths;
    std::string out;
    std::string model;
    std::string data;

    svm_options(int argc,char **argv)
    {
        putils::ProgramOptions opts;
        opts.addOptionsFromFile("svm.options");
        opts.parseCommandLine(argc,argv);
        if ( opts.hasValue("config")) {
            std::string config_file = opts.getValue("config");
            opts.parseFile(config_file.c_str());
        }
        out = opts.getValue("out");
        model = opts.getValue("model");
        data = opts.getValue("data");
        std::string task_str = opts.getValue("task");
        if ( task_str.compare("training")==0 ) {
            task = 0;
            eps = std::stod(opts.getValue("eps"));
            cost = std::stod(opts.getValue("cost"));
            kc1 = std::stod(opts.getValue("kc1"));
            kc2 = std::stod(opts.getValue("kc2"));
            kpow = std::stoi(opts.getValue("kpow"));
            ktype = std::stoi(opts.getValue("ktype"));
            kscale = std::stoi(opts.getValue("scale_kernel"));
            cache_mem = std::stoi(opts.getValue("cache_mem"));
            cache_size = std::stoi(opts.getValue("cache_size"));
            cache_precompute = std::stoi(opts.getValue("cache_precompute"));
            max_its = std::stoi(opts.getValue("max_its"));        
        }
        else
        {
            if ( task_str.compare("classify") == 0) {
                task_str = 1;
            }
            else
            {
                if ( task_str.compare("validate") == 0) {
                    task_str = 2;
                }
                else
                {
                    if ( task_str.compare("translate") == 0) {
                        task_str = 3;
                    }
                    else
                    {
                        std::cerr << "unknown task! " << task_str << "\n";
                        exit(-1);
                    }
                }
            }
        }
        nths = std::stoi(opts.getValue("nths"));
        write(std::cout);
    }
    std::ostream& write(std::ostream& os) const noexcept {
        os <<"svm options\n";
        os <<"data file  = "<<data<<"\n";
        os <<"model file = "<<model<<"\n";
        os <<"output     = "<<out<<"\n";
        os <<"# threads  = "<< nths << "\n";
        os <<"task       = "<< task << "\n";
        if ( task == 0) {
            os <<"scale kernel= "<< kscale << "\n";
            os <<"kernel type = "<< ktype << "\n";
            os <<"kernel pow  = "<<kpow << "\n";
            os <<"kernel c1   = "<<kc1 << "\n";
            os <<"kernel c2   = "<<kc2 << "\n";
            os <<"cache size  = "<<cache_size << "\n";
            os <<"cache mem   = "<<cache_mem << "\n";
            os <<"cache precomp = "<< cache_precompute << "\n";
            os <<"eps         = "<<eps << "\n";
            os <<"cost        = "<<cost << "\n";
            os <<"maxits      = "<<max_its << "\n";
        }
        return os;
    }

    friend std::ostream& operator << (std::ostream& os, const svm_options& opts)
    {
       return opts.write(os); 
    }
};
}
#endif