
#include "svm_classify.hpp"
#define MAX_NTHDS 64

void svm_classify(svm_options& options)
{
    stopwatch ctimer;
    ctimer.clear();
    ctimer.start();
    const double *vp;
    int i;
    double fx;
    svm_data data(options.data);
    const double *vecs = data.vecs;
    const double *y = data.y;
    int ntp=0;
    int nfp=0;
    int ntn=0;
    int nfn=0;
    int nth=0;
#ifdef _OPENMP
    int tid,xsz,sz0,sz,j;
    int fp,fn,tp,tn;
    int narr[MAX_NTHDS][4];
    int nth_i[1];

#pragma omp parallel private(tid,nth,xsz,sz,j,i,vp,fx,fp,fn,tp,tn) shared(narr,nth_i,vecs,y)
    {
        svm_model model(options.model);
        tid = omp_get_thread_num();
        nth = omp_get_num_threads();
        if (tid==0) nth_i[0]=nth;
        std::string my_file = options.out + "." + std::to_string(tid);
        std::ofstream out(my_file.c_str());
        xsz = data.nvecs % nth;
        sz = data.nvecs / nth;
        if (tid < xsz) {
            sz += 1;
            j = sz*tid;
        } else {
            j = sz*tid + xsz;
        }
        tp = 0;
        tn = 0;
        fp = 0;
        fn = 0;
        for (i=0; i<sz; ++i,++j) {
            vp = vecs + j * data.nfeat;
            fx = model.gensum(vecs+j*data.nfeat);
            if (fx > 0.) {
                if (y[j]>= 0.) {
                    ++tp;
                }else{
                    ++fp;
                }
            } else {
                if (y[j]>= 0.) {
                    ++tn;
                }else{
                    ++fn;
                }
            }
            out.write((char*)&y[j],sizeof(double));
            out.write((char*)&fx,sizeof(double));
        }
        out.close();
        narr[tid][0]=tp;
        narr[tid][1]=tn;
        narr[tid][2]=fp;
        narr[tid][3]=fn;
    }
    int nsum = 0;
    nth = nth_i[0];
    ntp = nfp = ntn = nfn = 0;
    for (i=0; i< nth; ++i) {
        ntp += narr[i][0];
        ntn += narr[i][1];
        nfp += narr[i][2];
        nfn += narr[i][3];
        nsum += narr[i][0] + narr[i][1] + narr[i][2] + narr[i][3];
    }
    std::cerr << "nth = " << nth << "\n";
    std::cerr << "nsum = " << nsum << "\n";
#else
    svm_model model(options.out);
    std::ofstream out;
    out.open(options.out.c_str());
    for (int i=0; i<data.nvecs; ++i) {
        fx = model.gensum(vecs+i*data.nfeat);
        if (fx >= 0.0) {
            if (y[i]>=0.0) ++ntp;
            else ++nfp;
        } else {
            if (y[i]>=0.0) ++ntn;
            else ++nfn;
        }
        out.write((char*)&y[i],sizeof(double));
        out.write((char*)&fx,sizeof(double));
    }
    out.close();
#endif
    ctimer.stop();
    std::cerr << "ntp = " << ntp << " nfp " << nfp << " ntn " << ntn << " nfn " << nfn << "\n";
    std::cerr << "classification time = " << ctimer.time() << " seconds\n";
    analyze(ntp,nfp,ntn,nfn);
}

