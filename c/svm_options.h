#ifndef SVM_OPTIONS_H
#define SVM_OPTIONS_H
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "utils.h"
#include "program_options.h"
#ifdef _OPENMP
#include <omp.h>
#endif

struct svm_options_t{
    double *vecs;
    double *y;
    double bias;
    double eps;
    double cost;
    double kc1;
    double kc2;
    int csize;
    int kpow;
    int ktype;
    int nvecs;
    int nfeat;
    int max_its;
    int task;
    int scale_kernel;
    int nths;
    char * out;
    char * model;
    char * data;
};

typedef struct svm_options_t svm_options_t;

svm_options_t * svm_options_init(int argc,char **argv);
void svm_options_free(svm_options_t* opts);
void svm_options_translate(svm_options_t* opts);
void svm_options_read_data_file(svm_options_t* opts);
void svm_options_write(svm_options_t* opts,FILE *fp);
#endif

