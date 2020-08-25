#ifndef SMO_SOLVER_H
#define SMO_SOLVER_H
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>
#include "stopwatch.h"
#include "svm_kernel_matrix.h"
#include "svm_options.h"
#include "utils.h"
#define USE_TIMERS

#ifdef _OPENMP
#define SVM_USE_OPENMP
#endif

typedef struct {
    svm_kernel_matrix_t *kmatrix;
    double **kmax;
    double **kmin;
    double *alfa;
    double *grad;
    double *y;
    int *status;
    double fun,bias,gap,cost,eps;
    int nvecs;
#ifdef USE_TIMERS
    stopwatch_t * kmat_timer;
    stopwatch_t * gap_timer;
    stopwatch_t * step_timer;
    stopwatch_t * grad_timer;
    stopwatch_t * find_timer;
#endif
} smo_solver_t;

struct  IndexPair {
    double value;
    int64_t index;
};

typedef struct IndexPair IndexPair;

IndexPair* PairMax(IndexPair *p1, IndexPair *p2);
IndexPair* PairMin(IndexPair *p1, IndexPair *p2);

void smo_solver_train(svm_options_t *opts);
int smo_solver_find_step(smo_solver_t *smo);
int smo_solver_take_step(smo_solver_t *smo,int imax,int imin);
void smo_solver_output_model_file(smo_solver_t *smo,svm_options_t *opts);
void smo_solver_find_gap(smo_solver_t *smo);
IndexPair * PairMax(IndexPair*,IndexPair*);
IndexPair * PairMin(IndexPair*,IndexPair*);
#endif