/*
 * SMOSolver.h
 *
 *  Created on: Jul 7, 2010
 *      Author: d3p708
 */
#include "smo_solver.h"

#define USE_TIMERS

void smo_solver_train(svm_options_t* options)
{
    int r,k,iter,maxits,nvecs;
    double diff=0.0;
    double fold=0.0;
    stopwatch_t *timer = timer_init();
    timer_clear(timer);
    smo_solver_t *smo = MALLOC_PTR(smo_solver_t);
    nvecs = smo->nvecs = options->nvecs;
    smo->alfa = (double*)Malloc(sizeof(double)*nvecs);
    smo->grad = (double*)Malloc(sizeof(double)*nvecs);
    smo->y = (double*)Malloc(sizeof(double)*nvecs);
    memcpy(smo->y,options->y,sizeof(double)*nvecs);
    maxits = options->max_its;
    smo->fun = 0.0;
    smo->gap = 0.0;
    smo->bias = 0.0;
    smo->cost = options->cost;
    smo->eps = options->eps;
    smo->status = (int*)malloc(sizeof(int)*nvecs);
    smo->kmax = (double**)Malloc(sizeof(double*));
    smo->kmin = (double**)Malloc(sizeof(double*));
    smo->kmatrix = svm_kernel_matrix_init(options);
#ifdef USE_TIMERS
    smo->gap_timer = timer_init();
    smo->step_timer = timer_init();
    smo->grad_timer = timer_init();
    smo->kmat_timer = timer_init();
    smo->find_timer = timer_init();
    timer_clear(smo->gap_timer);
    timer_clear(smo->step_timer);
    timer_clear(smo->grad_timer);
    timer_clear(smo->kmat_timer);
    timer_clear(smo->find_timer);
#endif
    memset(smo->alfa,0x0,sizeof(double)*nvecs);
    for ( k=0; k<nvecs; ++k) {
        smo->status[k] = -1;
        smo->y[k] = ( smo->y[k] > 0. ) ? 1.:-1.;
        smo->grad[k] = smo->y[k];
    }

    fprintf(stderr,"intialized smo solver %ld\n",maxits);
    timer_start(timer);

    for ( iter = 0; iter < maxits; ++iter ) {
        for ( k = 0; k < nvecs; ++k ) {
            r = smo_solver_find_step(smo);
            if ( r ) {
                fprintf(stderr," conv %d %d\n",iter,k);
                break;
            }
        }
        smo_solver_find_gap(smo);
        diff = smo->fun - fold;
        fold = smo->fun;
        fprintf(stderr, "Cycle          = %ld\n", iter );
        fprintf(stderr, "obj. function  = %20.10lf\n", smo->fun );
        fprintf(stderr, "diff. obj fun  = %le\n", diff );
        fprintf(stderr, "gap            = %le\n", smo->gap );
        fprintf(stderr, "bias           = %le\n\n", smo->bias);
        if ( r ) {
            fprintf(stderr, " converged! no more feasible step\n");
            break;
        }
        if ( smo->gap < smo->eps ) {
            fprintf(stderr, " converged! gap is within tolerance\n");
            break;
        }
    }
    timer_stop(timer);
    fprintf(stderr, "training time             = %le seconds\n",timer->elapsed_time);
#ifdef USE_TIMERS
    fprintf(stderr, "step time                 = %le seconds\n",smo->step_timer->elapsed_time);
    fprintf(stderr, "update time               = %le seconds\n",smo->find_timer->elapsed_time);
    fprintf(stderr, "gap time                  = %le seconds\n",smo->gap_timer->elapsed_time);
    fprintf(stderr, "grad time                 = %le seconds\n",smo->grad_timer->elapsed_time);
    fprintf(stderr, "kmat time                 = %le seconds\n",smo->kmat_timer->elapsed_time);
#endif
    smo_solver_output_model_file(smo,options);
}


void smo_solver_find_gap(smo_solver_t *smo) {
    int k;
    int nvecs = smo->nvecs;
    const double *alpha = smo->alfa;
    const double *grad = smo->grad;
    const double *dy = smo->y;
    const double zero = 0.0;
    double asum = 0.;
    double csum = 0.;
    int nfree=0;
    double fsum = 0.;
    double bsum = 0.;
#ifdef USE_TIMERS
    timer_start(smo->gap_timer);
#endif

#pragma omp parellel for  private(k) reduction(+,asum) reduction(+,fsum)
    for ( k = 0; k < nvecs; ++k ) {
        asum += alpha[k];
        fsum += alpha[k] * grad[k] * dy[k];
    }
    fsum = ( fsum + asum ) *0.5;
#pragma omp parellel for  private(k) reduction(+,bias) reduction(+,nfree)
    for ( k = 0; k < nvecs; ++k ) {
        if ( smo->status[k] != 0 )
            continue;
        bsum += grad[k] ;
        ++nfree;
    }
    bsum = -bsum;
///// if nfree = 0 then bias = 0 so no need for an else here
    if ( nfree )
        bsum /= nfree;
#pragma omp parellel for  private(k) reduction(+,csum)
    for ( k = 0; k < nvecs; ++k ) {
        csum += fmax( 0, ( dy[k] * ( grad[k] + bsum ) ) );
    }
    csum *= smo->cost;
    smo->gap = ( csum + asum - fsum - fsum ) / ( 1.0 + asum + csum - fsum );
    smo->fun = fsum;
    smo->bias = bsum;
#ifdef USE_TIMERS
    timer_stop(smo->gap_timer);
#endif
};

int smo_solver_take_step ( smo_solver_t *smo, int imax, int imin ) {
    double cost = smo->cost;
    const double *y = smo->y;
    double *grad = smo->grad;
    double ds, a1, a2, L, H, ai, aj, qc;
    int st1, st2, s, k;
    double da1,da2,gam;
    double zero = 0.;
    int nvecs = smo->nvecs;
#ifdef USE_TIMERS
    timer_start(smo->kmat_timer);
#endif
/*
    if (imax < 0 || imax>= smo->nvecs) {
        fprintf(stderr,"imax out of range %d\n",imax);
    }
    if (imin < 0 || imin>= smo->nvecs) {
        fprintf(stderr,"imin out of range %d\n",imin);
    }
*/
    svm_kernel_matrix_get_row(smo->kmatrix,imin,smo->kmin);
    svm_kernel_matrix_get_row(smo->kmatrix,imax,smo->kmax);
    const double *qmax= *(smo->kmax);
    const double *qmin= *(smo->kmin);
#ifdef USE_TIMERS
    timer_stop(smo->kmat_timer);
    timer_start(smo->step_timer);
#endif
    ds = y[imax] * y[imin];
    s = rint(ds);
    a1 = ai = smo->alfa[imax];
    a2 = aj = smo->alfa[imin];
    gam = a2 + ds * a1;
    if ( s > 0 ) {
        L = fmax(0.0,(gam-cost));
        H = fmin(cost,gam);
    } else {
        L = fmax(0.0,gam);
        H = fmin(cost,(gam+cost));
    }
    qc =fmax((qmax[imax] + qmin[imin] - qmax[imin] - qmax[imin]), TAU );
    a2 += ( ( grad[imin] - grad[imax] ) * y[imin] ) / qc;
    a2 = fmin( fmax ( a2, L ), H );
    a1 += ds * ( aj - a2 );
    if ( a1 > TAU ) {
        if ( a1 < ( cost - TAU ) ) {
            st1 = 0x0;
        } else {
            a1 = cost;
            a2 = gam - ds * cost;
            st1 = 0x1;
        }
    } else {
        a1 = zero;
        a2 = gam;
        st1 = -0x1;
    }
    if ( a2 > TAU ) {
        if ( a2 < ( cost - TAU ) ) {
            st2 = 0x0;
        } else {
            st2 = 0x1;
        }
    } else {
        st2 = -0x1;
    }
#ifdef USE_TIMERS
    timer_stop(smo->step_timer);
#endif

    if ( fabs ( a2 - aj ) > TAU ) {
        smo->alfa[imax] = a1;
        smo->alfa[imin] = a2;
        smo->status[imax] = st1;
        smo->status[imin] = st2;
#ifdef USE_TIMERS
        timer_start(smo->grad_timer);
#endif
        {
            da1 = ( y[imax] ) * ( a1 - ai );
            da2 = ( y[imin] ) * ( a2 - aj );
#pragma omp parellel for  private(k)  schedule(static,500)
            for ( k = 0; k < nvecs; ++k ) {
                grad[k] -= ( da1 * qmax[k] + da2 * qmin[k] );
            }
        }
#ifdef USE_TIMERS
        timer_stop(smo->grad_timer);
#endif
        return 0;
    }
    return -1;
};


IndexPair* PairMax(IndexPair *p1, IndexPair *p2) {
    if (p1->value >= p2->value) return p1;
    p1->value = p2->value;
    p1->index = p2->index;
    return p1;
}

IndexPair* PairMin(IndexPair *p1, IndexPair *p2) {
    if (p1->value <= p2->value) return p1;
    p1->value = p2->value;
    p1->index = p2->index;
    return p1;
}



#ifdef SVM_USE_OPENMP    
/*
int smo_solver_find_step(smo_solver_t *smo) {
    int k,i;
    double ys,gx;
    int nvecs = smo->nvecs;
    double max0 = 1.e300;
    IndexPair the_max;
    IndexPair the_min;
    IndexPair *my_max;
    IndexPair *my_min;
#ifdef USE_TIMERS
    timer_start(smo->find_timer);
#endif

#pragma omp declare reduction(MaxPair:IndexPair*:omp_out=PairMax(omp_out,omp_in))

#pragma omp declare reduction(MinPair:IndexPair*:omp_out=PairMin(omp_out,omp_in))

#pragma omp parallel 
  {
    the_max = &my_max;
    the_min = &my_min;
    the_max->value = -max0;
    the_max->index = -1;
    the_min->value = max0;
    the_min->index = -1;
#pragma omp parallel for private(k,gx,ys) reduction(PairMax:the_max) reduction(PairMin:the_min)
    for ( k = 0; k < nvecs; ++k ) {
        ys = smo->y[k] * smo->status[k];
        gx = smo->grad[k];
        if ( ys <= 0.0 ) {
            if (the_max->value < gx) {
                the_max->value = gx;
                the_max->index = k;
            }
        }
        if ( ys >= 0.0 ) {
            if (the_min->value > gx) {
                the_min->value = gx;
                the_min->index = k;
            }
        }
    }
  }
#ifdef USE_TIMERS
    timer_stop(smo->find_timer);
#endif
    if ( the_max->index != -1 || the_min->index != -1 ||
            the_max->index != the_min->index || fabs( the_max->value - the_min->value ) > smo->eps ) {
        return smo_solver_take_step (smo, the_max->index, the_min->index );
    }
    return -1;
}
*/

int smo_solver_find_step(smo_solver_t *smo) {
    int tid,nth;
    int64_t k,i;
    double ys,gx;
    double max0 = 1.e100;
    int nvecs = smo->nvecs;
    IndexPair the_max;
    IndexPair the_min;
    IndexPair max_pair[128];
    IndexPair min_pair[128];
    IndexPair my_max;
    IndexPair my_min;
#ifdef USE_TIMERS
    timer_start(smo->find_timer);
#endif    
#pragma omp parallel private(my_max,my_min,tid) shared(nth,max_pair,min_pair)
    {
        tid = omp_get_thread_num();
        nth = omp_get_num_threads();
        my_max.value = -max0;
        my_max.index = -1;
        my_min.value = max0;
        my_min.index = -1;
#pragma omp parallel for private(k,gx,ys)  
        for (k=0;k<nvecs;++k) {
            ys = smo->y[k] * smo->status[k];
            gx = smo->grad[k];
            if (ys <= 0.5) {
                if (my_max.value < gx) {
                    my_max.value = gx;
                    my_max.index = k;
                }
            } 
            if (ys >= -0.5) {
                if (my_min.value > gx) {
                    my_min.value = gx;
                    my_min.index = k;
                }
            }
        } 
        max_pair[tid].value = my_max.value;
        max_pair[tid].index = my_max.index;
        min_pair[tid].value = my_min.value;
        min_pair[tid].index = my_min.index;
    }
    the_max.value = max_pair[0].value;
    the_max.index = max_pair[0].index;
    the_min.value = min_pair[0].value;
    the_min.index = min_pair[0].index;
    for (i=1; i<nth; ++i) {
        if (the_max.value < max_pair[i].value) {
            the_max.value = max_pair[i].value;
            the_max.index = max_pair[i].index;
        }
        if (the_min.value < min_pair[i].value) {
            the_min.value = min_pair[i].value;
            the_min.index = min_pair[i].index;
        }
    }    
#ifdef USE_TIMERS
    timer_stop(smo->find_timer);
#endif
    if ( the_max.index != -1 || the_min.index != -1 ||
            the_max.index != the_min.index || fabs( the_max.value - the_min.value ) > smo->eps ) {
        return smo_solver_take_step (smo, (int)the_max.index, (int)the_min.index );
    }
    return -1;
}

#else
/*
 * Serial find step subroutine
 */
int smo_solver_find_step(smo_solver_t *smo) {
    int k,i;
    double ys,gx;
    int nvecs = smo->nvecs;
    double max0 = 1.e300;
    IndexPair the_max;
    IndexPair the_min;
    
#ifdef USE_TIMERS
    timer_start(smo->find_timer);
#endif
    the_max.value = -max0;
    the_max.index = -1;
    the_min.value = max0;
    the_min.index = -1;
    for ( k = 0; k < nvecs; ++k ) {
        ys = smo->y[k] * smo->status[k];
        gx = smo->grad[k];
        if ( ys <= 0. ) {
            if (the_max.value < gx) {
                the_max.value = gx;
                the_max.index = k;
            }
        }
        if ( ys >= 0. ) {
            if (the_min.value > gx) {
                the_min.value = gx;
                the_min.index = k;
            }
        }
    }
#ifdef USE_TIMERS
    timer_stop(smo->find_timer);
#endif
    if ( the_max.index != -1 || the_min.index != -1 ||
            the_max.index != the_min.index || fabs( the_max.value - the_min.value ) > smo->eps ) {
        return smo_solver_take_step (smo, the_max.index, the_min.index );
    }
    return -1;
}
#endif    


void smo_solver_output_model_file(smo_solver_t *smo, svm_options_t *opts) {
    const int nvecs = smo->nvecs;
    const int nfeat = opts->nfeat;
    int nsv = 0;
    int nbnd = 0;
    int ntp = 0 ;
    int nfp = 0 ;
    int ntn = 0 ;
    int nfn = 0 ;
    int itmp,k;
    double dtmp;
    FILE *out;
    FILE *out2;
    double fx;
    const double *alpha  = smo->alfa;
    const double *scal = smo->kmatrix->eval->scale;
    const double *y = smo->y;
    double *v;
#pragma omp parellel for  private(k) reduction(+,nsv) reduction(+,nbnd)
    for ( k = 0; k < nvecs; ++k ) {
        if ( smo->status[k] >= 0 ) {
            ++nsv;
            if ( smo->status[k] > 0 ) ++nbnd;
        }
    }
    fprintf(stderr, "# training vectors     = %d\n", smo->nvecs );
    fprintf(stderr, "# support vectors      = %d\n", nsv );
    fprintf(stderr, "# bound support vectors= %d\n", nbnd );
    out = Fopen(opts->model,"w");
    fwrite(&nsv,1,sizeof(int),out);
    fwrite(&(opts->nfeat),1,sizeof(int),out);
    fwrite(&(opts->ktype),1,sizeof(int),out);
    fwrite(&(opts->kpow),1,sizeof(int),out);
    fwrite(&(opts->kc1),1,sizeof(double),out);
    fwrite(&(opts->kc2),1,sizeof(double),out);
    fwrite(&(smo->bias),1,sizeof(double),out);
    fwrite(&(opts->scale_kernel),1,sizeof(int),out);
    for ( k = 0; k < nvecs; ++k ) {
        if ( smo->status[k] < 0 ) continue;
        dtmp = y[k] * scal[k] * alpha[k];
        fwrite(&dtmp,1,sizeof(double),out);
    }
    for ( k = 0; k < nvecs ; ++k) {
        if (smo->status[k]<0) continue;
        v = opts->vecs + k * nfeat;
        fwrite(v,nfeat,sizeof(double),out);
    }
    fclose(out);
    fprintf(stderr,"wrote model file\n");
    out2 = Fopen(opts->out,"w");
    for ( k = 0; k < nvecs; ++k ) {
        fx = ( smo->y[k] - smo->grad[k] - smo->bias );
        if ( fx > 0 ) {
            if ( smo->y[k] > 0. ) ++ntp;
            else ++nfp;
        } else {
            if ( smo->y[k] > 0. ) ++ntn;
            else ++nfn;
        }
        dtmp = smo->y[k];
        fwrite(&dtmp,1,sizeof(double),out2);
        fwrite(&fx,1,sizeof(double),out2);
    }
    fclose(out2);
    analyze ( ntp, nfp, ntn, nfn );
}
