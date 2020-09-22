#ifndef SVM_KERNEL_MATRIX_H
#define SVM_KERNEL_MATRIX_H

#include "svm_options.h"

typedef struct {
    const double * vecs;
    double * scale;
    double c1,c2;
    int kpow,ktype;
    int nvecs,nfeat;
} svm_kernel_eval_t;

typedef struct {
    svm_kernel_eval_t *eval;
    double *cache_row;
    int *cache_index;
    int nsize;
    int csize;
    int nvecs;
    int last;
} svm_kernel_matrix_t;

svm_kernel_eval_t * svm_kernel_eval_init(svm_options_t* opts);
void svm_kernel_eval_free(svm_kernel_eval_t * eval);
void svm_kernel_eval0(svm_kernel_eval_t *eval,double *row,int irow);
void svm_kernel_eval1(svm_kernel_eval_t *eval,double *row,int irow);
void svm_kernel_eval2(svm_kernel_eval_t *eval,double *row,int irow);
void svm_kernel_eval3(svm_kernel_eval_t *eval,double *row,int irow);
void svm_kernel_eval(svm_kernel_eval_t *eval,double *row,int irow);
svm_kernel_matrix_t * svm_kernel_matrix_init(svm_options_t* opts);
void svm_kernel_matrix_get_rows(svm_kernel_matrix_t *kmat,int imax,int imin,double **rmax,double **rmin);
void svm_kernel_matrix_free(svm_kernel_matrix_t *kmat);

#endif
