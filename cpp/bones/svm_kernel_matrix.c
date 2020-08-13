#include "svm_kernel_matrix.h"

inline svm_kernel_eval_t * svm_kernel_eval_init(const svm_options_t* opts)
{
    int i;
    double *scal;
    const double *vecs;
    double v,c1,c2;
    int nvecs,nfeat;

    svm_kernel_eval_t * eval = MALLOC_PTR(svm_kernel_eval_t);
    eval->vecs = opts->vecs;
    eval->nvecs = opts->nvecs;
    eval->nfeat = opts->nfeat;
    eval->ktype = opts->ktype;
    eval->kpow = opts->kpow;
    eval->c1 = opts->kc1;
    eval->c2 = opts->kc2;
    eval->scale = (double*)Malloc(sizeof(double)*(eval->nvecs));

    vecs = eval->vecs;
    c1 = eval->c1;
    c2 = eval->c2;
    nvecs = eval->nvecs;
    nfeat = eval->nfeat;

    if (opts->scale_kernel) {
        switch (eval->ktype) {
        case 0:
            for (i=0; i<nvecs; ++i) {
                v = eval0(vecs+i*nfeat,vecs+i*nfeat,nfeat);
                eval->scale[i] = 1./sqrt(v);
            }
            break;
        case 1:
            for (i=0; i<nvecs; ++i) {
                v = c1*eval0(vecs+i*nfeat,vecs+i*nfeat,nfeat)+c2;
                v = dpowi(v,eval->kpow);
                eval->scale[i] = 1./sqrt(v);
            }
            break;
        case 2:
            for (i=0; i<nvecs; ++i) {
                eval->scale[i] = 1.;
            }
            break;
        case 3:
            for (i=0; i<nvecs; ++i) {
                v = tanh( c1*eval0(vecs+i*nfeat,vecs+i*nfeat,nfeat)+c2 );
                eval->scale[i] = 1./sqrt(v);
            }
            break;
        default:
            fprintf(stderr,"unknown kernel type %d\n",eval->ktype);
            exit(EXIT_FAILURE);
        }
    } else {
        for (i=0; i<nvecs; ++i) eval->scale[i]=1.0;
    }

    fprintf(stderr,"Initialized Kernel Evaluator\n");
    return eval;
}

inline void svm_kernel_eval_free(svm_kernel_eval_t * eval)
{
    free(eval->scale);
    free(eval);
}

inline void svm_kernel_eval(svm_kernel_eval_t *eval,double *row,int irow)
{
    int i,j;
    int nvecs = eval->nvecs;
    int nfeat = eval->nfeat;
    const double *vecs = eval->vecs;
    const double *v1 = eval->vecs + irow*nfeat;
    const double *v2;
    double *s = eval->scale;
    double s1 = (eval->scale)[irow];
    double c1 = eval->c1;
    double c2 = eval->c2;
    int kpow = eval->kpow;
    double d,t;

    switch (eval->ktype) {
    case 0:
#pragma omp parallel for private(i)
        for (i=0; i<nvecs; ++i) row[i]=s1*s[i]*eval0(v1,vecs+i*nfeat,nfeat);
        return;
    case 1:
#pragma omp parallel for private(i)
        for (i=0; i<nvecs; ++i) row[i]=s1*s[i]*dpowi(c1*eval0(v1,vecs+i*nfeat,nfeat)+c2,kpow);
        return;
    case 2:
#pragma omp parallel for private(i,d,t,v2)     
        for (i=0; i<nvecs; ++i) {
            v2 = vecs + i * nfeat;
            d = 0.0;
            for (j=0; j<nfeat; ++j) {
                t = v1[j]-v2[j];
                d += t*t;
            }
            d *= c1;
            row[i]=exp(-d);
        }
        return;
    case 3:
#pragma omp parallel for private(i)
        for (i=0; i<nvecs; ++i) row[i]=s1*s[i]*tanh(c1*eval0(v1,vecs+i*nfeat,nfeat)+c2);
        return;
    }
}

inline svm_kernel_matrix_t * svm_kernel_matrix_init(svm_options_t* opts)
{
    int cache_size;
    int nvecs,i;

    svm_kernel_matrix_t *kmat = MALLOC_PTR(svm_kernel_matrix_t);

    kmat->eval = svm_kernel_eval_init(opts);
    kmat->csize = opts->csize;
    cache_size = kmat->csize;
    kmat->nvecs = opts->nvecs;
    nvecs = kmat->nvecs;
    kmat->last = 0;

    if (cache_size) {
        if (cache_size < 0 ) {
            fprintf(stderr,"cache size < 0\n");
            cache_size = nvecs/6;
            kmat->csize = cache_size;
        }
        kmat->cache_row = (double*)Malloc(sizeof(double)*nvecs*cache_size);
        kmat->cache_index = (int*)Malloc(sizeof(int)*cache_size);
        for (i=0; i<cache_size; ++i) (kmat->cache_index)[i]=-1;
    } else {
        kmat->cache_row = (double*)Malloc(sizeof(double)*nvecs*nvecs);
        kmat->cache_index = (int*)Malloc(sizeof(int)*nvecs);
#pragma omp parallel for private(i) 
        for (i=0; i<nvecs; ++i) {
            svm_kernel_eval(kmat->eval,kmat->cache_row+i*nvecs,i);
            kmat->cache_index[i] = i;
        }
        kmat->nsize = kmat->csize;
    }
    fprintf(stderr,"inited kmatrix # vecs = %d  cache size = %d \n",kmat->nvecs,kmat->csize);
    return kmat;
}

inline void svm_kernel_matrix_get_row(svm_kernel_matrix_t *kmat,int irow,double **row)
{
    int i,last,ix;
    int *kindex = kmat->cache_index;
    double *krows = kmat->cache_row;
    int nvecs = kmat->nvecs;
    int nsz= kmat->nsize;
    int cache_size = kmat->csize;
    if (kmat->csize) {
        ix = -1;
#pragma omp parallel for private(i) reduction(max:ix)
        for (i=0; i< kmat->csize; ++i) {
            if (kindex[i]==irow) {
                ix = i;
            }
        }
        if (ix>=0) {
            *row = krows + ix * nvecs;
            return;
        }
        // if we get here we didn't find row in cache
        last = kmat->last;
        svm_kernel_eval(kmat->eval,krows+last*nvecs,irow);
        *row = krows+last*nvecs;
        kmat->cache_index[last] = irow;
        kmat->last = (kmat->last + 1) % kmat->csize;
    } else {
        *row = krows+irow*nvecs;
    }
}

inline void svm_kernel_matrix_free(svm_kernel_matrix_t *kmat)
{
    free(kmat->eval);
    free(kmat->cache_row);
    free(kmat->cache_index);
    free(kmat);
}
