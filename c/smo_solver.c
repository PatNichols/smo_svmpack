/*
 * smo_solver.c
 *
 *  Created on: 2020
 *      Author: Patick Nichols
 */
#include "smo_solver.h"

smo_solver_t * smo_solver_init(svm_options_t *options) {
    int i;
    int nvecs = options->nvecs;
    smo_solver_t* smo = MALLOC_PTR(smo_solver_t);
    
    smo->nvecs = options->nvecs;
    smo->eps = options->eps;
    smo->cost = options->cost;
    smo->y = (double*)Malloc(sizeof(double)*options->nvecs);
    smo->alpha = (double*)Malloc(sizeof(double)*options->nvecs);
    smo->grad = (double*)Malloc(sizeof(double)*options->nvecs);
    smo->status = (int*)Malloc(sizeof(int)*options->nvecs);
    smo->fun = 0.0;
    smo->bias = 0.0;
    smo->gap = 1.e300;
    smo->index_max = (index_pair_t*)Malloc(sizeof(index_pair_t));
    smo->index_min = (index_pair_t*)Malloc(sizeof(index_pair_t));
    smo->kmax = (double**)Malloc(sizeof(double*));
    smo->kmin = (double**)Malloc(sizeof(double*));
    smo->kmatrix = svm_kernel_matrix_init(options);
    memset(smo->alpha,0x0,sizeof(double)*nvecs);
    for (i=0;i<nvecs;++i) {
        (smo->y)[i]=((options->y)[i] >= 0.) ? 1.0:-1.0;
        (smo->grad)[i] = (smo->y)[i];
        (smo->status)[i] = -1;
    }
#ifdef USE_TIMERS
    smo->gap_timer = stopwatch_init();
    smo->find_timer = stopwatch_init();
    smo->step_timer = stopwatch_init();
    smo->grad_timer = stopwatch_init();
    smo->kmat_timer = stopwatch_init();
#endif    
    return smo;
}

void smo_solver_free(smo_solver_t *smo)
{
#ifdef USE_TIMERS
    stopwatch_free(smo->kmat_timer);
    stopwatch_free(smo->grad_timer);
    stopwatch_free(smo->step_timer);
    stopwatch_free(smo->find_timer);
    stopwatch_free(smo->gap_timer);
#endif    
    svm_kernel_matrix_free(smo->kmatrix);
    free(smo->kmin);
    free(smo->kmax);
    free(smo->index_min);
    free(smo->index_max);
    free(smo->status);
    free(smo->grad);
    free(smo->alpha);
    free(smo->y);
    free(smo);
}

void smo_solver_train(svm_options_t * options)
{
    int r,k,iter;
    double diff=0.0;
    double fold=0.0;
    stopwatch_t * timer = stopwatch_init();
    smo_solver_t *smo = smo_solver_init(options);
    int nvecs = options->nvecs;
    int maxits = options->max_its;

    fprintf(stderr,"intialized smo solver %ld\n",maxits);
    stopwatch_start(timer);

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
        fprintf(stderr, "diff. obj fun  = %20.10le\n", diff );
        fprintf(stderr, "gap            = %20.10le\n", smo->gap );
        fprintf(stderr, "bias           = %20.10le\n\n", smo->bias);
        if ( r ) {
            fprintf(stderr, " converged! no more feasible step\n");
            break;
        }
        if ( smo->gap < smo->eps ) {
            fprintf(stderr, " converged! gap is within tolerance\n");
            break;
        }
    }
    stopwatch_stop(timer);
    fprintf(stderr, "training time             = %le seconds\n",timer->elapsed_time);
#ifdef USE_TIMERS
    fprintf(stderr, "step time                 = %le seconds\n",smo->step_timer->elapsed_time);
    fprintf(stderr, "update time               = %le seconds\n",smo->find_timer->elapsed_time);
    fprintf(stderr, "gap time                  = %le seconds\n",smo->gap_timer->elapsed_time);
    fprintf(stderr, "grad time                 = %le seconds\n",smo->grad_timer->elapsed_time);
    fprintf(stderr, "kmat time                 = %le seconds\n",smo->kmat_timer->elapsed_time);
#endif
    smo_solver_output_model_file(smo,options);
    smo_solver_free(smo);
}


void smo_solver_find_gap(smo_solver_t *smo) {
    int k;
    int nvecs = smo->nvecs;
    const double zero = 0.0;
    double asum = 0.;
    double csum = 0.;
    int nfree=0;
    double fsum = 0.;
    double bsum = 0.;
    const double *alfa = smo->alpha;
    const double *grad = smo->grad;
    const int *status = smo->status;
    const double *y = smo->y;
#ifdef USE_TIMERS
    stopwatch_start(smo->gap_timer);
#endif

#pragma omp parallel for  private(k) reduction(+:asum) reduction(+:fsum)
    for ( k = 0; k < nvecs; ++k ) {
        asum += alfa[k];
        fsum += alfa[k] * grad[k] * y[k];
    }
    fsum = ( fsum + asum ) * 0.5;
#pragma omp parallel for  private(k) reduction(+:bsum) reduction(+:nfree)
    for ( k = 0; k < nvecs; ++k ) {
        if ( status[k] == 0 ) {
            bsum += grad[k] ;
            ++nfree;
        }
    }
    bsum = -bsum;
///// if nfree = 0 then bias = 0 so no need for an else here
    if ( nfree )
        bsum /= nfree;
#pragma omp parallel for  private(k) reduction(+:csum)
    for ( k = 0; k < nvecs; ++k ) {
        csum += fmax( 0., ( y[k] * ( grad[k] + bsum ) ) );
    }
    csum *= smo->cost;
    smo->gap = ( csum + asum - fsum - fsum ) / ( 1.0 + asum + csum - fsum );
    smo->fun = fsum;
    smo->bias = bsum;
    fprintf(stderr,"asum = %lg csum = %lg\n",asum,csum);
#ifdef USE_TIMERS
    stopwatch_stop(smo->gap_timer);
#endif
};

inline int smo_solver_take_step (smo_solver_t *smo, int imax, int imin ) {
    double ds, a1, a2, L, H, ai, aj, qc;
    int st1, st2, s, k;
    double da1,da2,gam;
    double zero = 0.;
    const double tau = 2.e-14;
    double *y = smo->y;
    double *alpha = smo->alpha;
    int *status = smo->status;
    double *grad = smo->grad;
    const double *qmax;
    const double *qmin;
    const double cost = smo->cost;
#ifdef USE_TIMERS
    stopwatch_start(smo->kmat_timer);
#endif
#ifdef SVM_PARANOID
    if (imax < 0 || imax>= smo->nvecs) {
        fprintf(stderr,"imax out of range %d\n",imax);
        exit(EXIT_FAILURE);
    }
    if (imin < 0 || imin>= smo->nvecs) {
        fprintf(stderr,"imin out of range %d\n",imin);
        exit(EXIT_FAILURE);
    }
#endif
    svm_kernel_matrix_get_rows(smo->kmatrix,imax,imin,smo->kmax,smo->kmin);
    qmax= *(smo->kmax);
    qmin= *(smo->kmin);
#ifdef USE_TIMERS
    stopwatch_stop(smo->kmat_timer);
    stopwatch_start(smo->step_timer);
#endif
    ds = y[imax] * y[imin];
    s = rint(ds);
    a1 = ai = alpha[imax];
    a2 = aj = alpha[imin];
    gam = a2 + ds * a1;
    if ( s > 0 ) {
        L = fmax(0.0,(gam-cost));
        H = fmin(cost,gam);
    } else {
        L = fmax(0.0,gam);
        H = fmin(cost,(gam+cost));
    }
    qc =fmax((qmax[imax] + qmin[imin] - qmax[imin] - qmax[imin]), tau );
    a2 += ( ( grad[imin] - grad[imax] ) * y[imin] ) / qc;
    a2 = fmin( fmax ( a2, L ), H );
    a1 += ds * ( aj - a2 );
    if ( a1 > tau ) {
        if ( a1 < ( cost - tau ) ) {
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
    if ( a2 > tau ) {
        if ( a2 < ( cost - tau ) ) {
            st2 = 0x0;
        } else {
            st2 = 0x1;
        }
    } else {
        st2 = -0x1;
    }
#ifdef USE_TIMERS
    stopwatch_stop(smo->step_timer);
#endif
    if ( fabs ( a2 - aj ) > tau ) {
        alpha[imax] = a1;
        alpha[imin] = a2;
        status[imax] = st1;
        status[imin] = st2;
#ifdef USE_TIMERS
        stopwatch_start(smo->grad_timer);
#endif
        {
            const int nvecs = smo->nvecs;
            const double da1 = ( y[imax] ) * ( a1 - ai );
            const double da2 = ( y[imin] ) * ( a2 - aj );
#pragma omp parallel for  private(k) firstprivate(da1,da2)
            for ( k = 0; k < nvecs; ++k ) {
                grad[k] -= ( da1 * qmax[k] + da2 * qmin[k] );
            }
        }
#ifdef USE_TIMERS
        stopwatch_stop(smo->grad_timer);
#endif
        return 0;
    }
    return -1;
};

index_pair_t * index_pair_max(index_pair_t *p1,index_pair_t *p2) {
    if (p1->index != -1 && p1->value >= p2->value) return p1;
    p1->index = p2->index;
    p1->value = p2->value;
    return p1;
}

index_pair_t * index_pair_min(index_pair_t *p1,index_pair_t *p2) {
    if (p1->index != -1 && p1->value <= p2->value) return p1;
    p1->index = p2->index;
    p1->value = p2->value;
    return p1;
}

index_pair_t * index_pair_max_value(index_pair_t *p1,double v,int i)
{
    if (p1->index != -1 && p1->value >= v) return p1;
    p1->index = i;
    p1->value = v;
    return p1;
}

index_pair_t * index_pair_min_value(index_pair_t *p1,double v,int i)
{
    if (p1->index != -1 && p1->value <= v) return p1;
    p1->index = i;
    p1->value = v;
    return p1;
}

index_pair_t * index_pair_init(double v,int i)
{
    index_pair_t * p1 = (index_pair_t*)Malloc(sizeof(index_pair_t));
    p1->index = i;
    p1->value = v;
    return p1;
}

void index_pair_free(index_pair_t *p) {
    free(p);
}

#ifdef _OPENMP
/*
int smo_solver_find_step(smo_solver_t *smo) {
    int k,i;
    double ys,gx;
    const double eps_ = 1.e-14;
    const double max0 = 1.e300;
    double *grad = smo->grad;
    double *y = smo->y;
    int *status = smo->status;
    double max_val[129];
    double min_val[129];
    double the_max;
    double the_min;
    int max_index[129];
    int min_index[129];
    int the_max_index;
    int the_min_index;
    int sz,xsz,nth,tid;   
#ifdef USE_TIMERS
    stopwatch_start(smo->find_timer);
#endif
#pragma omp parallel private(tid,nth,sz,xsz,k,i,ys) shared(max_val,min_val,max_index,min_index)
    {
        tid = omp_get_thread_num();
        nth = omp_get_num_threads();
        sz = nvecs/nth;
        xsz = nvecs%nth;
        if (tid < xsz) {
            ++sz;
            k = tid*sz;
        }else{
            k = tid*sz + xsz;
        }
        max_val[tid]= -max0;
        max_index[tid] = -1;
        min_val[tid]= max0;
        min_index[tid] = -1;
        for (i = 0; i < sz; ++i) {
            ys = y[k] * status[k];
            if ( ys <= 0.  && max_val[tid] < grad[k]) {
                  max_val[tid] = grad[k];
                  max_index[tid] = k;
            }
            if ( ys >= 0. && min_val[tid] > grad[k] ) {
                  min_val[tid] = grad[k];
                  min_index[tid] = k;
            } 
            ++k;
        }
#pragam omp single 
        {
            the_max = max_val[0];
            the_max_index = max_index[0];
            for (i=0;i<nth;++i) {
                if (max_index[i]!=-1 && the_max < max_val[i]) {
                    the_max = max_val[i];
                    the_max_index = max_index[i];
                }
            }
        }
#pragma omp single
        {        
            the_min = min_val[0];
            the_min_index = min_index[0];
            for (i=0;i<nth;++i) {
                if (min_index[i]!=-1 && the_min > min_val[i]) {
                    the_min = min_val[i];
                    the_min_index = min_index[i];
                }
            }
        }
    }
#ifdef USE_TIMERS
    stopwatch_stop(smo->find_timer);
#endif
    if ( the_max_index != -1 || the_min_index != -1 ||
            the_max_index != the_min_index || fabs( the_max - the_min ) > eps_ ) {
        return smo_solver_take_step (smo,the_max_index, the_min_index );
    }
    return -1;
}
*/
int smo_solver_find_step(smo_solver_t *smo) {
    int k,i;
    double ys,gx;
    const double eps_ = 1.e-14;
    const double max0 = 1.e300;
    double *grad = smo->grad;
    double *y = smo->y;
    int *status = smo->status;
    index_pair_t * the_max = smo->index_max;
    index_pair_t * the_min = smo->index_min;
#ifdef USE_TIMERS
    stopwatch_start(smo->find_timer);
#endif
    the_max->value = -max0;
    the_max->index = -1;
    the_min->value = max0;
    the_min->index = -1;

#pragma omp declare reduction(MaxPair:index_pair_t*:omp_out=index_pair_max(omp_out,omp_in)) initializer(omp_priv=index_pair_init(-1.e300,-1))

#pragma omp declare reduction(MinPair:index_pair_t*:omp_out=index_pair_min(omp_out,omp_in)) initializer(omp_priv=index_pair_init(1.e300,-1))

#pragma omp parallel for private(k,gx,ys) reduction(MaxPair:the_max) reduction(MinPair:the_min)
    for ( k = 0; k < smo->nvecs; ++k ) {
        ys = y[k] * status[k];
        if ( ys <= 0. ) {
            the_max = index_pair_max_value(the_max,grad[k],k);
        }
        if ( ys >= 0. ) {
            the_min = index_pair_min_value(the_min,grad[k],k);
        }
    }
#ifdef USE_TIMERS
    stopwatch_stop(smo->find_timer);
#endif
    if ( the_max->index != -1 || the_min->index != -1 ||
            the_max->index != the_min->index || fabs( the_max->value - the_min->value ) > eps_ ) {
        return smo_solver_take_step (smo,the_max->index, the_min->index );
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
    double eps = 1.e-12;
    double max0 = 1.e300;
    double *grad = smo->grad;
    double *y = smo->y;
    int *status = smo->status;
    double max_value;
    double min_value;
    int max_index;
    int min_index;
#ifdef USE_TIMERS
    stopwatch_start(smo->find_timer);
#endif
    max_value = -max0;
    max_index = -1;
    min_value = max0;
    min_index = -1;

#ifdef USE_TIMERS
    stopwatch_start(smo->find_timer);
#endif
    for ( k = 0; k < smo->nvecs; ++k ) {
        ys = y[k] * status[k];
        gx = grad[k];
        if ( ys <= 0. ) {
            if (max_value < gx) {
                max_value = gx;
                max_index = k;
            }
        }
        if ( ys >= 0. ) {
            if (min_value > gx) {
                min_value = gx;
                min_index = k;
            }
        }
    }
#ifdef USE_TIMERS
    stopwatch_stop(smo->find_timer);
#endif
    if ( max_index != -1 || min_index != -1 ||
            max_index != min_index || fabs( max_value - min_value ) > eps ) {
        return smo_solver_take_step (smo,max_index,min_index);
    }
    return -1;
}
#endif


void smo_solver_output_model_file(smo_solver_t *smo,svm_options_t* opts) {
    const int nfeat = opts->nfeat;
    const int nvecs = smo->nvecs;
    int itmp;
    int nsv = 0;
    int nbnd = 0;
    int ntp = 0 ;
    int nfp = 0 ;
    int ntn = 0 ;
    int nfn = 0 ;
    int ntsv = 0;
    int nfsv = 0;
    double dtmp;
    double fx;
    const double *scal = smo->kmatrix->eval->scale;
    double *v;
    FILE *out;
    FILE *fout;
    int k;
    
#pragma omp parallel for  private(k) reduction(+:nsv) reduction(+:nbnd)
    for ( k = 0; k < nvecs; ++k ) {
        if ( smo->status[k] >= 0 ) {
            ++nsv;
            if ( smo->status[k] > 0 ) ++nbnd;
            if ( smo->y[k] >= 0) ++ntsv;
            else ++nfsv;
        }
    }

    fprintf(stderr,"# training vectors = %d\n",nvecs);
    fprintf(stderr,"# suppoer vectors  = %d\n",nsv);
    fprintf(stderr,"# bound support vectors = %d\n",nbnd);
    fprintf(stderr,"# true support vectors  = %d\n",ntsv);
    fprintf(stderr,"# false support vectors = %d\n",nfsv);

    out = Fopen(opts->model,"w");
    
    Fwrite((void *)&nsv,sizeof(int),1,out);
    Fwrite((void *)&nfeat,sizeof(int),1,out);
    Fwrite((void *)&(opts->ktype),sizeof(int),1,out);
    Fwrite((void *)&(opts->kpow),sizeof(int),1,out);
    Fwrite((void *)&(opts->kc1),sizeof(double),1,out);
    Fwrite((void *)&(opts->kc2),sizeof(double),1,out);
    Fwrite((void *)&(smo->bias),sizeof(double),1,out);
    itmp = 0;
    if (opts->scale_kernel) itmp=1;
    Fwrite((void *)&itmp,sizeof(int),1,out);
    for ( k = 0; k < nvecs; ++k ) {
        if (smo->status[k]<0) continue;
        dtmp = smo->y[k] * scal[k] * smo->alpha[k];
        Fwrite((void*)&dtmp,sizeof(double),1,out);
    }
    for ( k = 0; k < nvecs ; ++k) {
        if (smo->status[k]<0) continue;
        v = opts->vecs + k * nfeat;
        Fwrite((void*)v,sizeof(double),nfeat,out);
    }
    fclose(out);
    fprintf(stderr,"wrote model file\n");
    fout = Fopen(opts->out,"w");
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
        Fwrite((void*)&dtmp,sizeof(double),1,fout);
        Fwrite((void*)&fx,sizeof(double),1,fout);
    }
    fclose(fout);
    analyze ( ntp, nfp, ntn, nfn );
}

