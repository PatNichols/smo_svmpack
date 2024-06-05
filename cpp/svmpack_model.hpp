#ifndef SVM_CLASSIFY_HPP
#define SVM_CLASSIFY_HPP
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
namespace svmpack {
struct svm_model {
    double *yalf;
    double *vecs;
    double c1,c2,bias;
    int nfeat,nvecs;
    int ktype,kpower,kscale,owns;

    svm_model() = delete;
    ~svm_model()
    {
        delete [] vecs;
        delete [] yalf;
    }

    svm_model(const std::string& mod)
    {
        std::ifstream in(mod.c_str());
        if (!in) {
            std::cerr << "Could not open the model file " << mod << "\n";
            exit(EXIT_FAILURE);
        }
        in.read((char*)&nvecs,sizeof(int));
        in.read((char*)&nfeat,sizeof(int));
        in.read((char*)&ktype,sizeof(int));
        in.read((char*)&kpower,sizeof(int));
        in.read((char*)&c1,sizeof(double));
        in.read((char*)&c2,sizeof(double));
        in.read((char*)&bias,sizeof(double));
        in.read((char*)&kscale,sizeof(int));
        yalf = new double[nvecs];
        vecs = new double[nvecs*nfeat];
        in.read((char*)yalf,sizeof(double)*nvecs);
        in.read((char*)vecs,sizeof(double)*nvecs*nfeat);
        in.close();
        if (ktype>3 || ktype<0) {
            std::cerr << "unknown kernel type " << ktype << "\n";
            exit(EXIT_FAILURE);
        }
        std::cerr << "svm model file read " << mod << "\n";
        std::cerr << "nvecs_model = " << nvecs << "\n";
        std::cerr << "nfeat_model = " << nfeat << "\n";
        int nt = 0;
        int nf = 0;
        for (int i=0; i<nvecs; ++i) {
            if (yalf[i]>=0) ++nt;
            else ++nf;
        }
        std::cerr << "# true = " << nt << " #false = " << nf << "\n";
        owns = true;
    }
    constexpr double eval_dot(const double *v1,const double *v2) const noexcept {
        double s{0.0};
        for (int i=0; i<nfeat; ++i) s+=v1[i]*v2[i];
        return s;
    };
    constexpr double eval_diff(const double *v1,const double *v2) const noexcept {
        double s{0.0};
        for (int i=0; i<nfeat; ++i) {
            double t =v1[i]-v2[i];
            s+=t*t;
        }
        return s;
    };
    double scale_factor(const double *vp) const noexcept {
        if (!kscale || ktype==2) return 1.0;
        double s;
        switch (ktype) {
        case 0:
            s = c1 * eval_dot(vp,vp) + c2;
            return s;
        case 1:
            s = std::pow(c1*eval_dot(vp,vp)+c2,kpower);
            return s;
        case 2:
            return 1.0;
        case 3:
            s = std::tanh(c1*eval_dot(vp,vp)+c2);
            return s;
        }
        if ( s > tau ) {
            s = 1./std::sqrt(s);
        }else{
            s = 1.0;
        }
        return s;
    }
    double score(const double *vp) const noexcept {
        int i;
        double scal_p = scale_factor(vp);
        double sum = 0.0;
        switch (ktype) {
        case 0:
            for (i=0; i<nvecs; ++i) {
                sum += (c1 * eval_dot(vp,vecs+i*nfeat) + c2) * yalf[i];
            }
            return sum*scal_p-bias;
        case 1:
            for (i=0; i<nvecs; ++i) {
                sum += std::pow(c1*eval_dot(vp,vecs+i*nfeat)+c2,kpower) * yalf[i];
            }
            return sum*scal_p-bias;
        case 2:
            for (i=0; i<nvecs; ++i) {
                sum += std::exp( -(c1 * eval_diff(vp,vecs+i*nfeat) + c2) )* yalf[i];
            }
            return sum-bias;
        case 3:
            for (i=0; i<nvecs; ++i) {
                sum += std::tanh(c1*eval_dot(vp,vecs+i*nfeat)+c2) * yalf[i];
            }
            return sum*scal_p-bias;
        }
        return -1.e100;
    }
};
} // end namespace
#endif