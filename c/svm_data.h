#ifndef SVM_DATA_H
#define SVM_DATA_H
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "utils.h"

typedef struct {
    double *vecs;
    double *y;
    int nvecs;
    int nfeat;
} svm_data;

typedef struct {
    double *vecs;
    double *y;
    double c1,c2,bias;
    int nvecs,nfeat,ktype,kpow,kscal;
} svm_model;

svm_data * svm_data_init(const char *data_file);
void svm_data_free(svm_data *data);

#define svm_data_nvecs(data) (data)->nvecs
#define svm_data_nfeat(data) (data)->nfeat
#define svm_data_y(data) data->y
#define svm_data_vecs(data) (data)->vecs

svm_model * svm_model_init(const char *model_file);
void svm_model_free(svm_model *model);
double svm_model_gensum(svm_model *model,const double *vp);

#define svm_model_nvecs(data) (data)->nvecs
#define svm_model_nfeat(data) (data)->nfeat
#define svm_model_y(data) data->y
#define svm_model_vecs(data) (data)->vecs

#endif