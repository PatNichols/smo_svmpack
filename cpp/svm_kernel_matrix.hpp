#ifndef SVM_KERNEL_MATRIX_HPP
#define SVM_KERNEL_MATRIX_HPP
#include "stopwatch.hpp"
#include "svm_options.hpp"

struct svm_kernel_eval {
    const double * vecs;
    double * scale;
    double c1,c2;
    int kpow,ktype;
    int nvecs,nfeat;

    svm_kernel_eval(const svm_options& opts):vecs(opts.vecs),scale(0x0),c1(opts.kc1),c2(opts.kc2),
        kpow(opts.kpow),ktype(opts.ktype),nvecs(opts.nvecs),nfeat(opts.nfeat)
    {
        int i,j;
        scale = new double[nvecs];
        double sum,t;
        const double *v1;
        const double eps = 1.e-15;

        if (ktype>3 || ktype<0) {
            std::cerr << "Unknown Kernel Type " << ktype << "\n";
            exit(EXIT_FAILURE);
        }

        if (ktype==2) c1 = -c1;
        if (opts.scale_kernel) {
            switch (ktype) {
            case 0:
#pragma omp parallel for private(i,j,sum,v1)
                for (i=0; i<nvecs; ++i) {
                    v1 = vecs + i * nfeat;
                    sum =  0.;
                    for (j=0; j<nfeat; ++j) sum+=v1[j]*v1[j];
                    if (sum>eps) {
                        scale[i] = 1.0/sqrt(sum);
                    } else {
                        scale[i] = 0.0;
                    }
                }
                break;
            case 1:
#pragma omp parallel for private(i,j,sum,v1)
                for (i=0; i<nvecs; ++i) {
                    v1 = vecs + i * nfeat;
                    sum =  0.;
                    for (j=0; j<nfeat; ++j) sum+=v1[j]*v1[j];
                    sum = c1 * sum + c2;
                    sum = dpowi(sum,kpow);
                    if (sum>eps) {
                        scale[i] = 1.0/sqrt(sum);
                    } else {
                        scale[i] = 0.0;
                    }
                }
                break;
            case 2:
#pragma omp parallel for private(i) schedule(static,512)
                for (i=0; i<nvecs; ++i) {
                    scale[i]=1.0;
                }
                break;
            case 3:
#pragma omp parallel for private(i,j,sum,v1)
                for (i=0; i<nvecs; ++i) {
                    v1 = vecs + i * nfeat;
                    sum =  0.;
                    for (j=0; j<nfeat; ++j) sum+=v1[j]*v1[j];
                    sum = c1*sum+c2;
                    sum = tanh(sum);
                    if (sum>eps) {
                        scale[i] = 1.0/sqrt(sum);
                    } else {
                        scale[i] = 0.0;
                    }
                }
                break;
            }

        } else {
            for (i=0; i<nvecs; ++i) scale[i]=1.;
        }
        std::cerr << "kernel eval inited ";
        std::cerr << "ktype = " << ktype;
        std::cerr << " nvecs = " << nvecs << " nfeat = " << nfeat << "\n";
    }

    ~svm_kernel_eval()
    {
        delete [] scale;
    }

    void eval(double *row,int irow) {
        int i,j;
        double sum,t;
        double s1= scale[irow];
        const double *v1;
        const double *v2;
        switch (ktype) {
        case 0:
#pragma omp parallel for private(i,j,sum,v2,v1)
            for (i=0; i<nvecs; ++i) {
                v1 = vecs + irow*nfeat;
                v2 = vecs + i * nfeat;
                sum =  0.;
                for (j=0; j<nfeat; ++j) sum+=v1[j]*v2[j];
                row[i]=s1*scale[i]*sum;
            }
            return;
        case 1:
#pragma omp parallel for private(i,j,sum,v2,v1)
            for (i=0; i<nvecs; ++i) {
                v1 = vecs + irow*nfeat;
                v2 = vecs + i * nfeat;
                sum =  0.;
                for (j=0; j<nfeat; ++j) sum+=v1[j]*v2[j];
                sum = c1 * sum + c2;
                sum = dpowi(sum,kpow);
                row[i]=s1*scale[i]*sum;
            }
            return;
        case 2:
#pragma omp parallel for private(i,t,sum,v2,v1) schedule(static,512)
            for (i=0; i<nvecs; ++i) {
                v1 = vecs + irow*nfeat;
                v2 = vecs + i*nfeat;
                sum =  0.;
                for (j=0; j<nfeat; ++j) {
                    t = v1[j] - v2[j];
                    sum+=t*t;
                }
                row[i] = exp( sum * c1 );
            }
            return;
        case 3:
#pragma omp parallel for private(i,j,sum,v1,v2)
            for (i=0; i<nvecs; ++i) {
                v1 = vecs + irow*nfeat;
                v2 = vecs + i * nfeat;
                sum =  0.;
                for (j=0; j<nfeat; ++j) sum+=v1[j]*v2[j];
                sum = c1*sum+c2;
                sum = tanh(sum);
                row[i]=s1*scale[i]*sum;
            }
            return;
        }
    }
};

struct svm_kernel_matrix {
    svm_kernel_eval keval;
    double *cache_rows;
    int *cache_index;
    int nsize;
    int csize;
    int nvecs;
    int last;

    svm_kernel_matrix(const svm_options& opts):keval(opts),
        cache_rows(),cache_index(),csize(opts.csize),
        nvecs(opts.nvecs),last(0)
    {
        int i;
        if (csize==-1) {
            csize = nvecs/6;
        }
        if (csize==0) {
            cache_rows = new double[nvecs*nvecs];
            stopwatch stimer;
            stimer.clear();
#pragma omp parallel for private(i) //schedule(static,512)
            for (i=0; i<nvecs; ++i) keval.eval(cache_rows+i*nvecs,i);
            stimer.stop();
            std::cerr << "time to compute whole kernel matrix is : " << stimer.time() << "\n";
        } else {
            cache_rows = new double[csize*nvecs];
            cache_index = new int[csize];
            for (i=0; i<csize; ++i) cache_index[i]=-1;
        }
        std::cerr << "kernel matrix initialized ";
        std::cerr << "cache_size = " << csize << " nvecs = " << nvecs << " nfeat =" << opts.nfeat << "\n";
    }
    ~svm_kernel_matrix() {
        delete [] cache_index;
        delete [] cache_rows;
    }
    void get_row(int irow,double **R) {
        if (csize==0) {
            *R = cache_rows+irow*nvecs;
            return;
        }
        int i;
        int ifnd = -1;
#pragma omp parallel for private(i) reduction(max:ifnd)
        for (i=0; i<csize; ++i) {
            if (irow==cache_index[i]) ifnd = i;
        }
        if (ifnd>=0) {
            *R = cache_rows+ifnd*nvecs;
        } else {
            *R = cache_rows + last * nvecs;
            keval.eval(*R,irow);
            cache_index[last] = irow;
            last = (last + 1) % csize;
        }
    }

    inline void get_row(int imax,int imin,double **rmax,double **rmin) {
        int i;

        if (csize == 0) {
            *rmax = cache_rows + imax * nvecs;
            *rmin = cache_rows + imin * nvecs;
            return;
        }
        int imax_fnd = -1;
        int imin_fnd = -1;
#pragma omp parallel for private(i) shared(cache_index) reduction(max:imax_fnd) reduction(max:imin_fnd) // schedule(static,512)
        for (i=0; i<csize; ++i) {
            if (imax==cache_index[i]) imax_fnd = i;
            if (imin==cache_index[i]) imin_fnd = i;
        }
        if (imax_fnd!=-1) {
            *rmax = cache_rows+imax_fnd*nvecs;
        } else {
            // the next line is needed to prevent an overwrite if last==imin_fnd 
            if (last == imin_fnd) last = (last + 1)%csize;
            imax_fnd = last;
            keval.eval(cache_rows+last*nvecs,imax);
            *rmax = cache_rows + last * nvecs;
            cache_index[last] = imax;
            last = (last + 1) % csize;
        }
        if (imin_fnd!=-1) {
            *rmin = cache_rows+imin_fnd*nvecs;
        } else {
            // the next line is needed to prevent an overwrite if last==imax_fnd 
            if (last==imax_fnd) last=(last+1)%csize;
            keval.eval(cache_rows+last*nvecs,imin);
            *rmin = cache_rows + last * nvecs;
            cache_index[last] = imin;
            last = (last + 1) % csize;
        }
    }
};

#endif
