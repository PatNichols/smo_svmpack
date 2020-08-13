

#include "svm_classify.h"

svm_kfun_eval_t *svm_kfun_eval_init(svm_options_t* opts)
{
    svm_kfun_eval_t * k = MALLOC_PTR(svm_kfun_eval_t);
    k->vecs = opts->vecs;
    k->scale = opts->y;
    k->nvecs = opts->nvecs;
    k->nfeat = opts->nfeat;
    k->ktype = opts->ktype;
    k->kpow = opts->kpow;
    k->c1 = opts->kc1;
    k->c2 = opts->kc2;
    if (k->ktype==2) k->c1 = - k->c1;
    return k;
}

inline double svm_kfun_gensum(svm_kfun_eval_t *kfun,
                              const double *vp, int kscale) {
    int k,i;
    int nvecs = kfun->nvecs;
    int nfeat = kfun->nfeat;
    int kt = kfun->ktype;
    int kpow = kfun->kpow;
    const double *vecs = kfun->vecs;
    const double *scal = kfun->scale;
    double c1 = kfun->c1;
    double c2 = kfun->c2;
    double sum=0.;
    double sx= 1.0;

    if (kt!=2 && kscale) {
        switch (kt) {
        case 0:
            sx = eval0(vp,vp,nfeat);
            sx = 1./sqrt(sx);
            break;
        case 1:
            sx = c1*eval0(vp,vp,nfeat)+c2;
            sx = dpowi(sx,kpow);
            sx = 1./sqrt(sx);
            break;
        case 2:
            sx=1.0;
            break;
        case 3:
            sx = c1*eval0(vp,vp,nfeat)+c2;
            sx = tanh(sx);
            sx = 1./sqrt(sx);
            break;
        }
    }
    switch (kt) {
    case 0:
        for (k=0; k<nvecs; ++k) {
            sum += scal[k] * eval0(vp,vecs+k*nfeat,nfeat);
        }
        return sx*sum;
    case 1:
        for (k=0; k<nvecs; ++k) {
            sum += scal[k] * dpowi(c1*eval0(vp,vecs+k*nfeat,nfeat)+c2,kpow);
        }
        return sx*sum;
    case 2:
        for (k=0; k<nvecs; ++k) {
            sum += scal[k] * exp(c1 * eval1(vp,vecs+k*nfeat,nfeat));
        }
        return sx*sum;
    case 3:
        for (k=0; k<nvecs; ++k) {
            sum += scal[k] * tanh(c1*eval0(vp,vecs+k*nfeat,nfeat)+c2);
        }
        return sx*sum;
    }
    return -1.0;
}

void svm_classify(svm_options_t *opts)
{
    FILE *fp;
    char *new_name;
    char *data = opts->data;
    int nvecs,nfeat,i;
    double sx;
    double *y;
    double *vecs;
    double fx;
    size_t p;
    int off,xsz0,nsz0,nsz,tid,nth;
    int ntp,nfp,ntn,nfn;
    double *vp;
    double  bias = opts->bias;
    int kscale = opts->scale_kernel;
    svm_kfun_eval_t * kfun = svm_kfun_eval_init(opts);

    if (strstr(data,".tdo")==0x0) {
        new_name = strcat(data,".tdo");
        // this is not a tdo file
        translate_to_tdo(data,new_name);
    } else {
        new_name = data;
    }
    fp = Fopen(new_name,"r");
    fread(&nvecs,1,sizeof(double),fp);
    fread(&nfeat,1,sizeof(double),fp);
    y = (double*)Malloc(sizeof(double)*nvecs);
    fread(y,nvecs,sizeof(double),fp);
    p = ftell(fp);
    ntp = 0;
    nfp = 0;
    ntn = 0;
    nfn = 0;
#ifdef SVM_USE_OPENMP
    #pragma omp parallel private(tid,nsz0,xsz0,nsz,off,p,vecs,i,fx) \
    reduction(+:ntp) reduction(+:nfp) reduction(+:nfn) reduction(+:ntn)
    {
        tid = omp_get_thread_num();
        nth = omp_get_num_threads();
        nsz0 = nvecs/nth;
        xsz0 = nvecs%nth;
        if (tid < xsz0) {
            off = tid * (nsz0 + 1);
            nsz = nsz0 + 1;
        } else {
            off = tid * nsz0 + xsz0;
            nsz = nsz0;
        }
        vecs = (double*)Malloc(sizeof(double)*nfeat*nsz);
        off *= nfeat*sizeof(double);
        off += p;
        nsz0= nsz * nfeat;
        fseek(fp,off,SEEK_SET);
        Fread(vecs,nsz*nfeat,sizeof(double),fp);
        for (i=0; i<nsz; ++i) {
            vp = vecs + i * nfeat;
            fx = svm_kfun_gensum(kfun,vp,kscale) - bias;
            if (fx > 0.) {
                if (y[i]>0.) ++ntp;
                else ++nfp;
            } else {
                if (y[i]>0.) ++ntn;
                else ++nfn;
            }
        }
        free(vecs);
    }
#else
    vecs = (double*)Malloc(sizeof(double)*nfeat*nvecs);
    Fread(vecs,nvecs*nfeat,sizeof(double),fp);
    for (i=0; i<nvecs; ++i) {
        vp = vecs + i * nfeat;
        fx = svm_kfun_gensum(kfun,vp,kscale) - bias;
        if (fx > 0.) {
            if (y[i]>0.) ++ntp;
            else ++nfp;
        } else {
            if (y[i]>0.) ++ntn;
            else ++nfn;
        }
    }
    free(vecs);
#endif
    analyze(ntp,nfp,ntn,nfn);
    free(y);
    fclose(fp);
    svm_kfun_eval_free(kfun);
}
