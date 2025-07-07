#ifndef SVM_MODEL_H
#define SVM_MODEL_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "svm_options.h"
#include "svm_utils.h"
#include "svm_data.h"
#include "stopwatch.h"
#include "svm_fun.h"

struct svm_model {
    double *vecs;
    double *yalf;
    double c1,c2;
    double bias;
    int ktype;
    int kpow;
    int kscal;
    int nvecs;
    int nfeat;
};

typedef struct svm_model svm_model;

void svm_model_free(svm_model *model);
svm_model * svm_model_init(const char * model_file,int tid);
double svm_model_scale_factor(svm_model* model,const double *vp);
double svm_model_kernel_sum(svm_model* model,const double *vp);
void svm_model_classify(svm_options_t *options);
#endif

