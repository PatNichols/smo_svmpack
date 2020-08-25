
#include "svm_classify.hpp"


void svm_classify(svm_options& options)
{
        stopwatch ctimer;
        ctimer.clear();
        ctimer.start();
        int i;
        double fx;
        svm_model model_i;
        svm_model model(options.model);
        svm_data data(options.data);
        int ntp=0;
        int nfp=0;
        int ntn=0;
        int nfn=0;
/*
#ifdef _OPENMP
        int tid,nth,xsz,sz0,sz,off,j;   
        int fp = 0;
        int fn = 0;
        int tn = 0;
        int tp = 0;
#pragma omp parallel private(tid,nth,xsz,sz0,sz,off,model_i,fp,fn,tp,tn) 
      {
        fp = fn = tp = tn = 0;
        model_i = model;
        tid = omp_get_thread_num(); 
        nth = omp_get_num_threads();
        std::string my_file = options.out + "." + std::to_string(tid);
        std::ofstream out(my_file.c_str());
        xsz = data.nvecs % nth;
        sz0 = data.nvecs / nth;
        if (tid < xsz) {
            sz = sz0 + 1;
            off = sz*tid;
        }else{
            sz = sz0;
            off = sz*tid + xsz;
        }
        if (tid==0) {
            std::cerr << "tid = " << tid << " nth " << nth << " nsz = " << sz << "\n";
            std::cerr << "file size = " << (sz*2*sizeof(double)) << "\n";
        }
#pragma omp parallel for private(i,j,model_i,fx)
        for (i=0;i<sz;++i) {
            j = off + i;
            fx = model_i.gensum(data.vecs+j*data.nfeat);
            if (fx >= 0.0) {
                if (data.y[j]>=0.0) ++tp;
                else ++fp;
            }else{
                if (data.y[j]>=0.0) ++tn;
                else ++fn;
            } 
            out.write((char*)&(data.y[j]),sizeof(double));
            out.write((char*)&fx,sizeof(double));
        }
        out.close();
#pragma omp critical
        {
            ntp += tp;
            nfp += fp;
            ntn += tn;
            nfn += fn;
        }
      }
#else 
*/
        std::ofstream out;
        out.open(options.out.c_str());
        for (int i=0;i<data.nvecs;++i) {
            fx = model.gensum(data.vecs+i*data.nfeat);
            if (fx >= 0.0) {
                if (data.y[i]>=0.0) ++ntp;
                else ++nfp;
            }else{
                if (data.y[i]>=0.0) ++ntn;
                else ++nfn;
            } 
            out.write((char*)&(data.y[i]),sizeof(double));
            out.write((char*)&fx,sizeof(double));
        }
        out.close();
//#endif
        ctimer.stop();
        std::cerr << "ntp = " << ntp << " nfp " << nfp << " ntn " << ntn << " nfn " << nfn << "\n";
        std::cerr << "classification time = " << ctimer.time() << " seconds\n";
        analyze(ntp,nfp,ntn,nfn);
}

