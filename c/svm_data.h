#ifndef SVM_DATA_H
#define SVM_DATA_H
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "svm_utils.h"

typedef struct {
    double *vecs;
    double *y;
    int nvecs;
    int nfeat;
} svm_data;

svm_data * svm_data_init(const char *data_file);
void svm_data_free(svm_data *data);

#define svm_data_nvecs(data) (data)->nvecs
#define svm_data_nfeat(data) (data)->nfeat
#define svm_data_y(data) data->y
#define svm_data_vecs(data) (data)->vecs
#endif