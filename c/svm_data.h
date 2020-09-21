#ifndef SVM_DATA_H
#define SVM_DATA_H

#include "utils.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef struct {
    double *vecs;
    double *y;
    int nvecs;
    int nfeat;
} svm_data;

typedef struct {
    double *vecs;
    double *yalfa;
    double c1,c2,bias;
    int nvecs,nfeat,ktype,kpow,kscal;
} svm_model;

svm_data * svm_data_init(const char *data_file);
void svm_data_free(svm_data *data);
void svm_data_read_libsvm(svm_data *data,const char *data_file);

svm_model * svm_model_init(const char *model_file);
void svm_model_free(svm_model *model);
double svm_model_gensum(svm_model *model,const double *vp);
#endif