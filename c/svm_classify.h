#ifndef SVM_CLASSIFY_H
#define SVM_CLASSIFY_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#define SVM_USE_OPENMP
#endif
#include "utils.h"
#include "svm_options.h"

typedef struct
{
    double *vecs;
    double *scale;
    double c1,c2;
    int nvecs,nfeat,ktype,kpow;
} svm_kfun_eval_t;

svm_kfun_eval_t *svm_kfun_eval_init(svm_options_t* opts);

#define svm_kfun_eval_free(k) free((k))

double svm_kfun_gensum(svm_kfun_eval_t *kfun,const double *vp,int scale_k);

void svm_classify(svm_options_t *opts);

#endif
