#ifndef SVM_OPTIONS_H
#define SVM_OPTIONS_H
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include "utils.hpp"
#include "program_options.hpp"

struct svm_options {
    double *vecs;
    double *y;
    double bias;
    double eps;
    double cost;
    double kc1;
    double kc2;
    int csize;
    int kpow;
    int ktype;
    int nvecs;
    int nfeat;
    int max_its;
    int task;
    int scale_kernel;
    int nths;
    std::string out;
    std::string model;
    std::string data;

    svm_options(int argc,char **argv);

    ~svm_options() {
        delete [] vecs;
        delete [] y;
    }

    void read_libsvm_data_file() noexcept;
    void read_tdo_data_file() noexcept;
    void read_model_file() noexcept;
    void read_data_file() noexcept;
    void write_tdo_file(std::string&);
    void write_libsvm_file(std::string&);
    void translate();
    void write_file(FILE *fp) const noexcept;

    std::ostream& write(std::ostream& os) const noexcept {
        double vsize = double(sizeof(double)* nvecs * nfeat)/1048576.0;
        double ksize = double(sizeof(double)* nvecs * csize)/1048576.0;

        os <<"svm options\n";
        os <<"data file = "<<data<<"\n";
        os <<"model file= "<<model<<"\n";
        os <<"output    = "<<out<<"\n";
        os <<"# vectors  = "<< nvecs << "\n";
        os <<"# features = "<< nfeat << "\n";
        os <<"# threads  = "<< nths << "\n";
        os <<"scale kernel= "<<scale_kernel << "\n";
        os <<"kernel type = "<< ktype << "\n";
        os <<"kernel pow  = "<<kpow << "\n";
        os <<"kernel c1   = "<<kc1 << "\n";
        os <<"kernel c2   = "<<kc2 << "\n";
        os <<"task        = "<<task << "\n";
        if (task==0) {
            os <<"cache size  = "<<csize << "\n";
            os <<"eps         = "<<eps << "\n";
            os <<"cost        = "<<cost << "\n";
            os <<"maxits      = "<<max_its << "\n";
        }
        if (task==1) {
            os <<"bias        = "<<bias << "\n";
        }
        os << "vecs size   = "<<vsize <<" MB\n";
        os << "kmat size   = "<<ksize <<" MB\n";
        return os;
    }
};

inline std::ostream& operator << (std::ostream& os,const svm_options& opts) {
    return opts.write(os);
}

#endif