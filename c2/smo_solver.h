#ifndef smo_solver_h
#define smo_solver_h
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include "svm_utils.h"
#include "svm_fun.h"
#include "svm_options.h"
#include "svm_kernel_matrix.h"
#include "stopwatch.h"

#define SVM_USE_TIMERS 1

struct index_pair_t {
    double value;
    int index;
};

typedef struct index_pair_t index_pair_t;

struct smo_solver_t {
    double *y;
    double *alpha;
    double *grad;
    int *status;
    double fun;
    double gap;
    double eps;
    double cost;
    double bias;
    index_pair_t *index_max;
    index_pair_t *index_min;
    svm_kernel_matrix_t * kmatrix;
    double **kmax;
    double **kmin;
    int nvecs;
    int maxits;
#ifdef SVM_USE_TIMERS
    stopwatch_t *kmat_timer;
    stopwatch_t *grad_timer;
    stopwatch_t *gap_timer;
    stopwatch_t *find_timer;
    stopwatch_t *step_timer;
#endif
};

typedef struct smo_solver_t smo_solver_t;

void smo_solver_train(svm_options_t * options);

index_pair_t * index_pair_max(index_pair_t* p1,index_pair_t *p2);
index_pair_t * index_pair_min(index_pair_t* p1,index_pair_t *p2);
index_pair_t * index_pair_max_value(index_pair_t *p,double v,int i);
index_pair_t * index_pair_min_value(index_pair_t *p,double v,int i);
index_pair_t * index_pair_init(double v,int i);
void index_pair_free(index_pair_t *p);

smo_solver_t * smo_solver_init(svm_options_t *options);
void smo_solver_free(smo_solver_t *smo);
void smo_solver_find_gap(smo_solver_t *smo);
int smo_solver_take_step (smo_solver_t *smo, int imax, int imin );
int smo_solver_find_step(smo_solver_t *smo);
void smo_solver_output_model_file(smo_solver_t *smo,svm_options_t* opts);
#endif
