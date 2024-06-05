#include "svm_classify.hpp"
namespace svmpack
{
void svm_validate(const svm_options& options)
{
    svm_model model(options.model);
    long ntp(0);
    long nfp(0);
    long ntn(0);
    long nfn(0);
#ifdef _OPENMP
// read input and score in parallel
    std::ifstream in(options.data);
    if (!in) {
        std::cerr << "could not open the input file " << options.data << "\n";
        exit(-1);
    }
    in.read((char*)&nvecs,sizeof(int));
    in.read((char*)&nfeat,sizeof(int));
    size_t off1 = sizeof(int);
    off1 = off1 + off1;
#pragma omp parallel reduction(+,ntp) reduction(+,nfp) reduction(+,ntn) reduction(+,nfn)
    {
        int tid = omp_thread_num();
        int nth = omp_num_threads();
        long sz0 = nvecs/nth;
        long xsz = nvecs%nth;
        // read in my labels 
        long len = sz0;
        long off = sz0 * tid;
        if (tid < xsz) {
            len += 1;
            off += tid;
        }else{
            off += xsz;
        }
        double * y = new double[len];
        double * vecs = new double[nfeat*len];
        size_t yoff = off1 + off * sizeof(double);
        in.seekg(yoff,ios::beg);
        in.read((char*)&y,sizeof(double)*len);
        size_t voff = off1 + sizeof(double)*(nvecs + off * nfeat);
        in.seekg((voff,ios::beg);        
        in.read((char*)vecs,len*nfeat*sizeof(double));
        for (auto k=0;k<len;++k)
        {
            double fx = model.score(vecs+k*nfeat);
            double yx = y[k];
            if ( fx > -1.e-16) {
                if ( yx > -1.e-16) {
                    ++ntp;
                }else{
                    ++nfp;
                }
            }else{
                if ( yx > -1.e-16) {
                    ++nfn;
                }else{
                    ++ntn;
                }
            }
        }    
        delete [] vecs;
        delete [] y;
    }
    in.close();
#else
    svm_data data(options.data);
    const double * vecs = data.vectors();
    const double * y = data.labels();
    int nvecs = data.num_vectors();
    int nfeat = data.num_features();
    int nfeat_mod = model.num_features();
    // if the # of features of the model > # of features of data
    // there will be a oob in score below
    if ( nfeat < nfeat_mod) {
        data.edit(nfeat_mod);
        nfeat = nfeat_mod;
    }
    for (auto k=0;k<nvecs;++k)
    {
        double fx = model.score(vecs+k*nfeat);
        double yx = y[k];
        if ( fx > -1.e-16) {
            if ( yx > -1.e-16) {
                ++ntp;
            }else{
                ++nfp;
            }
        }else{
            if ( yx > -1.e-16) {
                ++nfn;
            }else{
                ++ntn;
            }        
        }
    }
#endif    
    analyze(ntp,ntn,nfp,nfn);
}

void svm_classify(const svm_options& options)
{
    std::cerr << "classification not implemented yet!\n";
}
}