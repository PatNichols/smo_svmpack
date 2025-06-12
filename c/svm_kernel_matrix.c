#include <math.h>
#include <omp.h>
#include "svm_kernel_matrix.h"
#include "svm_fun.h"

double kern_powi(double x,int m) {
    double w,y,z;
    if (m<0) {
        m = -m;
        x = 1./x;
    }
    switch(m) {
    case 0:
        return 1.0;
    case 1:
        return x;
    case 2:
        return x*x;
    case 3:
        return x*x*x;
    case 4:
        y =x*x;
        return y*y;
    case 5:
        y = x*x;
        return y*y*x;
    case 6:
        y = x*x;
        return y*y*y;
    case 7:
        y = x*x;
        return y*y*y*x;
    case 8:
        y= x*x;
        y= y*y;
        return y*y;
    case 9:
        y = x*x;
        y = y*y;
        return y*y*x;
    case 10:
        y = x*x;
        z = y*y;
        z = y*y;
        return z*y;
    case 11:
        y = x*x;
        z = y*y;
        z = y*y;
        return z*y*x;
    case 12:
        y = x*x;
        z = y*y*y;
        return z*z;
    default:
        y = x*x;
        z = y*y*y;
        z = z*z;
        return z*kern_powi(x,(m-12));
    }
    return 1.e300;
}

svm_kernel_eval_t * svm_kernel_eval_init(const svm_options_t * opts) {
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

//    if (eval->ktype==2) eval->c1 = -1.* eval->c1;
    if (opts->scale_kernel && opts->ktype!=2) {
        switch (eval->ktype) {
        case 0:
            for (i=0; i<nvecs; ++i) {
                v = svm_dot(vecs+i*nfeat,vecs+i*nfeat,nfeat);
                eval->scale[i] = 1./sqrt(v);
            }
            break;
        case 1:
            for (i=0; i<nvecs; ++i) {
                v = c1*svm_dot(vecs+i*nfeat,vecs+i*nfeat,nfeat)+c2;
                v = kern_powi(v,eval->kpow);
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
                v = tanh( c1*svm_dot(vecs+i*nfeat,vecs+i*nfeat,nfeat)+c2 );
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

void svm_kernel_eval_free(svm_kernel_eval_t * eval) {
    free(eval->scale);
    free(eval);
}

void svm_kernel_eval(svm_kernel_eval_t *eval,double *row,int irow) {
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
        #pragma omp parallel for private(i,s1)
        for (i=0; i<nvecs; ++i) row[i]=s1*s[i]*svm_dot(v1,vecs+i*nfeat,nfeat);
        return;
    case 1:
        #pragma omp parallel for private(i,s1)
        for (i=0; i<nvecs; ++i) row[i]=s1*s[i]*kern_powi(c1*svm_dot(v1,vecs+i*nfeat,nfeat)+c2,kpow);
        return;
    case 2:
        #pragma omp parallel for private(i,j,d,t,v2) shared(row,vecs) schedule(static)
        for (i=0; i<nvecs; ++i) {
            v2 = vecs + i * nfeat;
            d = 0.0;
            for (j=0; j<nfeat; ++j) {
                t = v1[j]-v2[j];
                d += t*t;
            }
            row[i]=exp(-c1*d);
        }
        return;
    case 3:
        #pragma omp parallel for private(i)
        for (i=0; i<nvecs; ++i) row[i]=s1*s[i]*tanh(c1*svm_dot(v1,vecs+i*nfeat,nfeat)+c2);
        return;
    }
}

svm_kernel_matrix_t * svm_kernel_matrix_init(svm_options_t* opts) {
    int cache_size;
    int nvecs,i;

    svm_kernel_matrix_t *kmat = MALLOC_PTR(svm_kernel_matrix_t);

    kmat->eval = svm_kernel_eval_init(opts);
    kmat->csize = opts->csize;
    cache_size = kmat->csize;
    kmat->nvecs = opts->nvecs;
    nvecs = kmat->nvecs;
    kmat->last = 0;

    if (cache_size>=0) {
        if (cache_size < 0 ) {
            fprintf(stderr,"cache size < 0\n");
            cache_size = nvecs/6;
            kmat->csize = cache_size;
        }
        kmat->cache_row = (double*)Malloc(sizeof(double)*nvecs*cache_size);
        kmat->cache_index = (int*)Malloc(sizeof(int)*cache_size);
        for (i=0; i<cache_size; ++i) (kmat->cache_index)[i]=-1;
    } else {
        kmat->csize = nvecs;
        kmat->cache_row = (double*)Malloc(sizeof(double)*nvecs*nvecs);
        kmat->cache_index = (int*)Malloc(sizeof(int)*nvecs);
        for (i=0; i<nvecs; ++i) {
            svm_kernel_eval(kmat->eval,kmat->cache_row+i*nvecs,i);
            kmat->cache_index[i] = i;
        }
        kmat->nsize = kmat->csize;
    }
    fprintf(stderr,"inited kmatrix # vecs = %d  cache size = %d \n",kmat->nvecs,kmat->csize);
    return kmat;
}

void svm_kernel_matrix_get_rows(svm_kernel_matrix_t *kmat,
                                int imax,int imin,double **rmax,double **rmin) {
    int i,imax_f,imin_f;
    int *cache_index = kmat->cache_index;
    double *krows = kmat->cache_row;
    int nvecs = kmat->nvecs;
    int last = kmat->last;
    int cache_size = kmat->csize;
    if (cache_size < nvecs) {
        imax_f = -1;
        imin_f = -1;
        #pragma omp parallel for private(i) reduction(max:imax_f) reduction(max:imin_f) schedule(static)
        for (i=0; i< cache_size; ++i) {
            if (cache_index[i]==imax) imax_f = i;
            if (cache_index[i]==imin) imin_f = i;
        }
        if (imax_f>=0) {
            *rmax = krows + imax_f * nvecs;
        } else {
            // if we get here we didn't find row in cache
            if (imin_f == last) last = (last+1)%cache_size;
            svm_kernel_eval(kmat->eval,krows+last*nvecs,imax);
            *rmax = krows+last*nvecs;
            cache_index[last] = imax;
            last = (last + 1) % cache_size;
        }
        if (imin_f>=0) {
            *rmin = krows + imin_f * nvecs;
        } else {
            // if we get here we didn't find row in cache
            if (imax_f == last) last = (last+1)%cache_size;
            svm_kernel_eval(kmat->eval,krows+last*nvecs,imin);
            *rmin = krows+last*nvecs;
            cache_index[last] = imin;
            last = (last + 1) % cache_size;
        }
        kmat->last = last;
        return;
    } else {
        *rmax = krows+imax*nvecs;
        *rmin = krows+imin*nvecs;
    }
}

void svm_kernel_matrix_free(svm_kernel_matrix_t *kmat) {
    free(kmat->eval);
    free(kmat->cache_row);
    free(kmat->cache_index);
    free(kmat);
}
