#ifndef SMO_SOLVER_HPP
#define SMO_SOLVER_HPP
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>
#include "stopwatch.hpp"
#include "svm_kernel_matrix.hpp"
#include "svm_options.hpp"
#include "utils.hpp"
#define USE_TIMERS

#ifdef _OPENMP
#define SVM_USE_OPENMP
#endif

struct smo_solver {
    svm_kernel_matrix kmatrix;
    double **kmax;
    double **kmin;
    double *alfa;
    double *grad;
    double *y;
    int *status;
    double fun,bias,gap,cost,eps;
    int nvecs;
#ifdef USE_TIMERS
    stopwatch kmat_timer;
    stopwatch gap_timer;
    stopwatch step_timer;
    stopwatch grad_timer;
    stopwatch find_timer;
#endif
    smo_solver(svm_options& options):kmatrix(options),kmax(0x0),kmin(0x0),
        alfa(0x0),grad(0x0),y(options.y),
        status(0x0),fun(0.),
        bias(0.),gap(0.),cost(options.cost),eps(options.eps),
        nvecs(options.nvecs)
#ifdef USE_TIMERS    
        ,kmat_timer(),gap_timer(),step_timer(),grad_timer(),find_timer()
#endif
    {
    int k;
    alfa = new double[nvecs];
    grad = new double[nvecs];
    status = new int[nvecs];
    kmax = new double*();
    kmin = new double*();
    for (k=0; k<nvecs; ++k) {
        alfa[k]=0.;
    }
    for (k=0;k<nvecs;++k) {
        status[k] = -1;
    }
    for (k=0;k<nvecs;++k) {
        y[k] = ( y[k] >= 0. ) ? 1.:-1.;
        grad[k] = y[k];
    }
    }
    ~smo_solver() {
        delete kmin;
        delete kmax;
        delete [] status;
        delete [] grad;
        delete [] alfa;
    }
    int find_step() noexcept;
    int take_step(int,int) noexcept;
    void find_gap() noexcept;
    void output_model_file(const svm_options& opts);

};

struct  IndexPair {
    double value;
    int64_t index;
    
    IndexPair():value(0),index(-1) {}

    IndexPair(double p,int64_t i):value(p),index(i) {}
    
    IndexPair& Max(IndexPair& p) {
        if (value < p.value) {
            value = p.value;
            index = p.index;
        }
        return *this;
    }
    IndexPair& Min(IndexPair& p) {
        if (value > p.value) {
            value = p.value;
            index = p.index;
        }
        return *this;
    }
    IndexPair& Max(const double pvalue,const int64_t pindex) {
        if (value < pvalue) {
            value = pvalue;
            index = pindex;
        }
        return *this;
    }
    IndexPair& Min(const double pvalue,const int64_t pindex) {
        if (value > pvalue) {
            value = pvalue;
            index = pindex;
        }
        return *this;
    }   
};

void smo_solver_train(svm_options& opts);

#endif



