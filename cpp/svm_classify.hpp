#ifndef SVM_CLASSIFY_HPP
#define SVM_CLASSIFY_HPP

#include <cstdio>
#include <cstdlib>
#include <cmath>
#ifdef SVM_USE_OPENMP
#include <omp.h>
#endif
#include "utils.hpp"
#include "svm_options.hpp"

struct svm_eval
{
    const double *vecs;
    const double *yalfa;
    double c1,c2;
    int nvecs,nfeat,ktype,kpow,kscal;

    svm_eval(svm_options& options):vecs(options.vecs),yalfa(options.y),c1(options.kc1),c2(options.kc2),
        nvecs(options.nvecs),nfeat(options.nfeat),ktype(options.ktype),
        kpow(options.kpow),kscal(options.scale_kernel)
    {
        if (ktype>3 || ktype <0) {
            std::cerr << "unknown kernel type "<< ktype<< "\n";
            exit(EXIT_FAILURE);
        }                
    }
    
    ~svm_eval() {
        vecs = nullptr;
        yalfa = nullptr;
    }


    constexpr double eval0(const double *v1,const double *v2) const noexcept {
        double s(0.);
        for (int i=0;i<nfeat;++i) s+=v1[i]*v2[i];
        return s;
    }
    constexpr double eval2(const double *v1,const double *v2) const noexcept {
        double s(0);
        for (int i=0;i<nfeat;++i) {
            auto t=v1[i]-v2[i];
            s+=t*t;
        }
        return s;
    }

    inline double gensum(const double *vp) {
    int k,i;
    double sum=0.;
    double sx= 1.0;

    if (kscal) {
        switch (ktype) {
        case 0:
            sx = eval0(vp,vp);
            sx = 1./sqrt(sx);
            break;
        case 1:
            sx = c1*eval0(vp,vp)+c2;
            sx = dpowi(sx,kpow);
            sx = 1./sqrt(sx);
            break;
        case 2:
            break;
        case 3:
            sx = c1*eval0(vp,vp)+c2;
            sx = tanh(sx);
            sx = 1./sqrt(sx);
            break;
        }
    }
    switch (ktype) {
    case 0:
        for (k=0; k<nvecs; ++k) {
            sum += yalfa[k] * eval0(vp,vecs+k*nfeat);
        }
        return sx*sum;
    case 1:
        for (k=0; k<nvecs; ++k) {
            sum += yalfa[k] * dpowi(c1*eval0(vp,vecs+k*nfeat)+c2,kpow);
        }
        return sx*sum;
    case 2:
        sum = 0.0;
        for (k=0; k<nvecs; ++k) {
            sum += yalfa[k] * exp( -c1 * eval2(vp,vecs+k*nfeat));
        }
        return sum;
    case 3:
        for (k=0; k<nvecs; ++k) {
            sum += yalfa[k] * tanh(c1*eval0(vp,vecs+k*nfeat)+c2);
        }
        return sx*sum;
    }
    std::cerr << "kernel type is bad " << ktype << "\n";
    exit(EXIT_FAILURE);
    }
};

void svm_classify(svm_options& opts);
#endif