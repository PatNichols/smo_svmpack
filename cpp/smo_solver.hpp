#ifndef SMO_SOLVER_HPP
#define SMO_SOLVER_HPP
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <fstream>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "stopwatch.hpp"
#include "svm_kernel_matrix.hpp"
#include "svm_options.hpp"
#include "putils_cxx.hpp"
#include "svmpack_math.hpp"
#define SVMPACK_USE_TIMERS

namespace svmpack {

struct smo_solver {
    svm_kernel_matrix kmatrix;
    double *alfa;
    double *grad;
    const double *y;
    int *status;
    double fun,bias,gap,cost,eps;
    int nvecs;
#ifdef SVMPACK_USE_TIMERS
    stopwatch kmat_timer;
    stopwatch gap_timer;
    stopwatch step_timer;
    stopwatch grad_timer;
    stopwatch find_timer;
#endif
    smo_solver(const svm_options& opts, const svm_data& data):
        kmatrix(opts,data),
        alfa(0x0),
        grad(0x0),
        y(data.labels()),
        status(0x0),
        fun(0.),
        bias(0.),
        gap(0.),
        cost(opts.cost),
        eps(opts.eps),
        nvecs(data.num_vectors())
#ifdef SVMPACK_USE_TIMERS
        ,kmat_timer(),
        gap_timer(),
        step_timer(),
        grad_timer(),
        find_timer()
#endif
    {
        int k;
        alfa = new double[nvecs];
        grad = new double[nvecs];
        status = new int[nvecs];
        for (k=0; k<nvecs; ++k) {
            alfa[k]=0.;
        }
        for (k=0; k<nvecs; ++k) {
            status[k] = -1;
        }
        for (k=0; k<nvecs; ++k) {
            grad[k] = ( y[k] > -1.e-16 ) ? 1.:-1.;
        }
    }
    ~smo_solver() {
        delete [] status;
        delete [] grad;
        delete [] alfa;
    }
    int find_step() noexcept;
    int take_step(int,int) noexcept;
    void find_gap() noexcept;
    void output_model_file(const svm_options&,const svm_data&);
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

void smo_solver_train(const svm_options& opts);
}
#endif
