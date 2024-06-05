#ifndef SVM_DATA_H
#define SVM_DATA_H
#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include "putils_cxx.hpp"

namespace svmpack {
struct svm_data {
    double *vecs;
    double *y;
    int nvecs;
    int nfeat;
    
    svm_data():
        vecs(nullptr),
        y(nullptr),
        nvecs(0),
        nfeat(0)
    {}

    svm_data(const std::string& data_):
        vecs(nullptr),
        y(nullptr),
        nvecs(0),
        nfeat(0)
    {
        read_data_file(data_);
    }

    ~svm_data() {
        delete [] vecs;
        delete [] y;
    }

    void read_libsvm_data_file(const std::string&);
    void read_tdo_data_file(const std::string&);
    void read_data_file(const std::string&);
    void write_tdo_data_file(const std::string&);
    void write_libsvm_data_file(const std::string&, const double&);
    void write_data_file(const std::string&);
    
    int num_vectors() const noexcept { return nvecs;}
    int num_features() const noexcept { return nfeat;}
    const double * vectors() const noexcept { return vecs;}
    const double * labels() const noexcept { return y;}

    void edit(std::size_t nfeat_new);
};

void svm_translate(const std::string& data_in);
} // end namespace
#endif