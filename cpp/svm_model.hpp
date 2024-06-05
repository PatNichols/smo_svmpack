#ifndef SVM_MODEL_HPP
#define SVM_MODEL_HPP
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include "putils_cxx.hpp"
#include "svmpack_math.hpp"
namespace svmpack {
class svm_model
{
    double * vecs;
    double * yalf;
    double bias;
    double c1;
    double c2;
    int scale_kernel;
    int ktype;
    int kpow;
    int nvecs;
    int nfeat;
public:
    svm_model(const std::string& model)
    {
        read_model_file(model);
    }

    ~svm_model()
    {
        delete [] vecs;
        delete [] yalf;
    }

    constexpr int num_features() const noexcept { return nfeat;}

    void read_model_file(const std::string& model)
    {
        /// training task
        std::ifstream in;
        in.open(model.c_str());
        if (!in) {
            throw putils::FileOpenError(__FUNCTION__,model.c_str(),"r");
        }
        in.read((char*)&nvecs,sizeof(int));
        in.read((char*)&nfeat,sizeof(int));
        in.read((char*)&ktype,sizeof(int));
        in.read((char*)&kpow,sizeof(int));
        in.read((char*)&c1,sizeof(double));
        in.read((char*)&c2,sizeof(double));
        in.read((char*)&bias,sizeof(double));
        in.read((char*)&(scale_kernel),sizeof(int));
        vecs = new double[nvecs*nfeat];
        yalf = new double[nvecs];
        in.read((char*)yalf,nvecs*sizeof(double));
        in.read((char*)vecs,nvecs*nfeat*sizeof(double));
        in.close();
        std::cerr << "# model vectors = " << nvecs << "\n";
        std::cerr << "# model features= " << nfeat << "\n";
        std::cerr << "model kernel function = " << ktype << "\n";
        std::cerr << "      parameters " << ktype << " " << c1 << " " << c2 << "\n";
        std::cerr << "bias = " << bias << "\n";
        std::cerr << "scale kernel = " << scale_kernel << "\n";
    }

    double scale_vector(const double *vin) const noexcept
    {
        double s = 0.0;
        if (scale_kernel) {
            switch(ktype)
            {
            case 0:
                s = c1 * dot(vin,vin,nfeat) + c2;
                break;
            case 1:
                s = std::pow(c1 * dot(vin,vin,nfeat) + c2,kpow);
                break;
            case 2:
                s = 1.0;
                return s;
            case 3:
                s = tanh(c1*dot(vin,vin,nfeat)+c2);
                break;
            }
            if (fabs(s) > 1.e-16) {
                s = 1./sqrt(s);
            } else {
                s = 1.;
            }
            return s;
        }
        return 1.0;
    }

    double score(const double *vin) const noexcept {

        double scale_in = scale_vector(vin);
        double sum{0.0};
        switch (ktype)
        {
        case 0:
            for (int i=0; i<nvecs; ++i)
            {
                const double * v = vecs + i * nfeat;
                sum += yalf[i] * ( c1 * dot(v,vin,nfeat) + c2);
            }
            break;
        case 1:
            for (int i=0; i<nvecs; ++i)
            {
                const double * v = vecs + i * nfeat;
                sum += yalf[i] * std::pow( (c1 * dot(v,vin,nfeat) + c2), kpow);
            }
            break;
        case 2:
            for (int i=0; i<nvecs; ++i)
            {
                const double * v = vecs + i * nfeat;
                sum += yalf[i] * exp( - (c1 * diff_nrm2(v,vin,nfeat) + c2) );
            }
            break;
        case 3:
            for (int i=0; i<nvecs; ++i)
            {
                const double * v = vecs + i * nfeat;           
                sum += yalf[i] * tanh( c1 * dot(v,vin,nfeat) + c2);
            }
            break;
        default:
            break;
        }
        return scale_in*sum-bias;
    }
}; // end class
} // end namespace
#endif
