
#include "svm_classify.hpp"

#ifdef _OPENMP
#include "omp.h"
#endif

void svm_classify(svm_options& opts)
{
    std::string data(opts.data);
    int nvecs,nfeat,i,j,tid,nth;
    int nsz,sz,xsz,off;
    double *y;
    double fx;
    size_t p;
    int ntp,nfp,ntn,nfn;
    double *vp;
    const double  bias = opts.bias;
    std::ifstream in;
    std::ofstream out;
    
    std::string new_name;
    if (data.find(".tdo")==std::string::npos) {
        new_name = opts.data + ".tdo";
        opts.translate();
    }else{
        new_name = opts.data;
    }
    in.open(new_name.c_str());
    in.read((char*)&nvecs,sizeof(int));
    in.read((char*)&nfeat,sizeof(int));
    y  = new double[nvecs];
    in.read((char*)y,nvecs*sizeof(double));
    p = in.tellg();
    std::cerr << "Reading Data file\n";
    std::cerr << "# vecs = " << nvecs << " # features= " << nfeat << "\n";
    ntp = 0;
    nfp = 0;
    ntn = 0;
    nfn = 0;        
#ifdef _OPENMP
    in.close();
    std::ifstream inx;
#pragma omp parallel private(inx,tid,nsz,xsz,sz,off,p,i,j,vp,fx) shared(y,nvecs,opts)\
    reduction(+:ntp) reduction(+:nfp) reduction(+:nfn) reduction(+:ntn)
    {
        tid = omp_get_thread_num();
        nth = omp_get_num_threads();
        nsz = nvecs/nth;
        xsz = nvecs%nth;
        if (tid < xsz) {
            sz = nsz+1;
            off = sz*tid;
        }else{
            sz = nsz;
            off = sz*tid + xsz;
        }
        p = sizeof(int)*2 + sizeof(double)*(nvecs + off * nfeat);
        svm_eval kfun(opts);
        inx.open(opts.data.c_str());
        inx.seekg(p);
        double *vecs = new double[sz*nfeat];
        inx.read((char*)vecs,sizeof(double)*nfeat*sz);
        inx.close();
        j = off;
        for (i=0; i<sz; ++i,++j) {
            vp = vecs + i * nfeat;
            fx = kfun.gensum(vp)- bias;
            if (fx > 0.) {
                if (y[j]>= 0.) ++ntp;
                else ++nfp;
            } else {
                if (y[j]>= 0.) ++ntn;
                else ++nfn;
            }
        }
        delete [] vecs;
    }
#else
    svm_eval kfun(opts);
    double *vecs = new double[nvecs*nfeat];
    in.read((char*)vecs,sizeof(double)*nfeat*nvecs);
    in.close();
    for (i=0; i<nvecs; ++i) {
        vp = vecs + i * nfeat;
        fx = kfun.gensum(vp)- bias;
        if (fx > 0.) {
            if (y[i]> 0.) ++ntp;
            else ++nfp;
        } else {
            if (y[i]> 0.) ++ntn;
            else ++nfn;
        }
        out.write((char*)&y[i],sizeof(double));
        out.write((char*)&fx,sizeof(double));
    }
    delete [] vecs;
#endif
    delete [] y;
    analyze(ntp,nfp,ntn,nfn);
}
