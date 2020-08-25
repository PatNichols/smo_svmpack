/*
 * SMOSolver.h
 *
 *  Created on: Jul 7, 2010
 *      Author: d3p708
 */
#include "smo_solver.hpp"

#define USE_TIMERS

void smo_solver_train(svm_options& options)
{
    int r,k,iter,maxits,nvecs;
    double diff=0.0;
    double fold=0.0;
    stopwatch timer;
    timer.clear();
    smo_solver smo(options);
    nvecs = options.nvecs;
    maxits = options.max_its;
#ifdef USE_TIMERS
    smo.gap_timer.clear();
    smo.step_timer.clear();
    smo.grad_timer.clear();
    smo.kmat_timer.clear();
    smo.find_timer.clear();
#endif

    fprintf(stderr,"intialized smo solver %ld\n",maxits);
    timer.start();

    for ( iter = 0; iter < maxits; ++iter ) {
        for ( k = 0; k < nvecs; ++k ) {
            r = smo.find_step();
            if ( r ) {
                fprintf(stderr," conv %d %d\n",iter,k);
                break;
            }
        }
        smo.find_gap();
        diff = smo.fun - fold;
        fold = smo.fun;
        fprintf(stderr, "Cycle          = %ld\n", iter );
        fprintf(stderr, "obj. function  = %20.10lf\n", smo.fun );
        fprintf(stderr, "diff. obj fun  = %20.10le\n", diff );
        fprintf(stderr, "gap            = %20.10le\n", smo.gap );
        fprintf(stderr, "bias           = %20.10le\n\n", smo.bias);
        if ( r ) {
            fprintf(stderr, " converged! no more feasible step\n");
            break;
        }
        if ( smo.gap < smo.eps ) {
            fprintf(stderr, " converged! gap is within tolerance\n");
            break;
        }
    }
    timer.stop();
    fprintf(stderr, "training time             = %le seconds\n",timer.elapsed_time);
#ifdef USE_TIMERS
    fprintf(stderr, "step time                 = %le seconds\n",smo.step_timer.elapsed_time);
    fprintf(stderr, "update time               = %le seconds\n",smo.find_timer.elapsed_time);
    fprintf(stderr, "gap time                  = %le seconds\n",smo.gap_timer.elapsed_time);
    fprintf(stderr, "grad time                 = %le seconds\n",smo.grad_timer.elapsed_time);
    fprintf(stderr, "kmat time                 = %le seconds\n",smo.kmat_timer.elapsed_time);
#endif
    smo.output_model_file(options);
}


void smo_solver::find_gap() noexcept {
    int k;
    const double zero = 0.0;
    double asum = 0.;
    double csum = 0.;
    int nfree=0;
    double fsum = 0.;
    double bsum = 0.;
#ifdef USE_TIMERS
    gap_timer.start();
#endif

#pragma omp parellel for  private(k) reduction(+,asum) reduction(+,fsum)
    for ( k = 0; k < nvecs; ++k ) {
        asum += alfa[k];
        fsum += alfa[k] * grad[k] * y[k];
    }
    fsum = ( fsum + asum ) / 2.0;
#pragma omp parellel for  private(k) reduction(+,bias) reduction(+,nfree)
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
#pragma omp parellel for  private(k) reduction(+,csum)
    for ( k = 0; k < nvecs; ++k ) {
        csum += fmax( 0, ( y[k] * ( grad[k] + bsum ) ) );
    }
    csum *= cost;
    gap = ( csum + asum - fsum - fsum ) / ( 1.0 + asum + csum - fsum );
    fun = fsum;
    bias = bsum;
    std::cerr << "asum = " << asum << " csum  =" << csum << "\n";
#ifdef USE_TIMERS
    gap_timer.stop();
#endif
};

int smo_solver::take_step ( int imax, int imin ) noexcept {
    double ds, a1, a2, L, H, ai, aj, qc;
    int st1, st2, s, k;
    double da1,da2,gam;
    double zero = 0.;
    const double tau = 2.e-14;
#ifdef USE_TIMERS
    kmat_timer.start();
#endif
#ifdef SVM_PARANOID
    if (imax < 0 || imax>= nvecs) {
        fprintf(stderr,"imax out of range %d\n",imax);
        exit(EXIT_FAILURE);
    }
    if (imin < 0 || imin>= nvecs) {
        fprintf(stderr,"imin out of range %d\n",imin);
        exit(EXIT_FAILURE);
    }
#endif    
    kmatrix.get_row(imax,imin,kmax,kmin);
    const double *qmax= *kmax;
    const double *qmin= *kmin;
#ifdef USE_TIMERS
    kmat_timer.stop();
    step_timer.start();
#endif
    ds = y[imax] * y[imin];
    s = rint(ds);
    a1 = ai = alfa[imax];
    a2 = aj = alfa[imin];
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
    step_timer.stop();
#endif

    if ( fabs ( a2 - aj ) > tau ) {
        alfa[imax] = a1;
        alfa[imin] = a2;
        status[imax] = st1;
        status[imin] = st2;
#ifdef USE_TIMERS
        grad_timer.start();
#endif
        {
            const double da1 = ( y[imax] ) * ( a1 - ai );
            const double da2 = ( y[imin] ) * ( a2 - aj );
#pragma omp parellel for  private(k,da1,da2)  schedule(static,1000)
            for ( k = 0; k < nvecs; ++k ) {
                grad[k] -= ( da1 * qmax[k] + da2 * qmin[k] );
            }
        }
#ifdef USE_TIMERS
        grad_timer.stop();
#endif
        return 0;
    }
    return -1;
};


IndexPair& PairMax(IndexPair &p1, IndexPair &p2) noexcept {
    if (p1.index!=-1 && p1.value >= p2.value) return p1;
    p1.value = p2.value;
    p1.index = p2.index;
    return p1;
}

IndexPair& PairMin(IndexPair &p1, IndexPair &p2) noexcept {
    if (p1.index!=-1 && p1.value <= p2.value) return p1;
    p1.value = p2.value;
    p1.index = p2.index;
    return p1;
}



#ifdef SVM_USE_OPENMP
int smo_solver::find_step() noexcept {
    int k,i;
    double ys,gx;
    double max0 = 1.e300;
    IndexPair the_max;
    IndexPair the_min;
#ifdef USE_TIMERS
    find_timer.start();
#endif
    the_max.value = -max0;
    the_max.index = -1;
    the_min.value = max0;
    the_min.index = -1;

#pragma omp declare reduction(MaxPair:IndexPair:omp_out.Max(omp_in)) initializer(omp_priv=IndexPair(-1.e300,-1))

#pragma omp declare reduction(MinPair:IndexPair:omp_out.Min(omp_in)) initializer(omp_priv=IndexPair(1.e300,-1))

#pragma omp parallel for private(k,gx,ys) reduction(MaxPair:the_max) reduction(MinPair:the_min) schedule(static,1000)
    for ( k = 0; k < nvecs; ++k ) {
        ys = y[k] * status[k];
        gx = grad[k];
        if ( ys <= 0. ) {
            the_max.Max(grad[k],k);
        }
        if ( ys >= 0. ) {
            the_min.Min(grad[k],k);
        }
    }
#ifdef USE_TIMERS
    find_timer.stop();
#endif
    if ( the_max.index != -1 || the_min.index != -1 ||
            the_max.index != the_min.index || fabs( the_max.value - the_min.value ) > eps ) {
        return take_step (the_max.index, the_min.index );
    }
    return -1;
}
#else
/*
 * Serial find step subroutine
 */
int smo_solver::find_step() noexcept {
    int k,i;
    double ys,gx;
    double max0 = 1.e300;
    IndexPair the_max;
    IndexPair the_min;

#ifdef USE_TIMERS
    find_timer.start();
#endif
    the_max.value = -max0;
    the_max.index = -1;
    the_min.value = max0;
    the_min.index = -1;
    for ( k = 0; k < nvecs; ++k ) {
        ys = y[k] * status[k];
        gx = grad[k];
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
    find_timer.stop();
#endif
    if ( the_max.index != -1 || the_min.index != -1 ||
            the_max.index != the_min.index || fabs( the_max.value - the_min.value ) > eps ) {
        return take_step ((int)the_max.index,(int)the_min.index );
    }
    return -1;
}
#endif


void smo_solver::output_model_file(const svm_options& opts) {
    const int nfeat = opts.nfeat;
    int nsv = 0;
    int nbnd = 0;
    int ntp = 0 ;
    int nfp = 0 ;
    int ntn = 0 ;
    int nfn = 0 ;
    int k;
    double dtmp;
    double fx;
    const double *scal = kmatrix.keval.scale;
    double *v;
#pragma omp parellel for  private(k) reduction(+,nsv) reduction(+,nbnd)
    for ( k = 0; k < nvecs; ++k ) {
        if ( status[k] >= 0 ) {
            ++nsv;
            if ( status[k] > 0 ) ++nbnd;
        }
    }

    std::cerr << "# training vectors      = " << nvecs << "\n";
    std::cerr << "# support vectors       = " << nsv << "\n";
    std::cerr << "# bound support vecs    = " << nbnd << "\n";

    std::ofstream out;
    out.open(opts.model.c_str());
    out.write((char*)&nsv,sizeof(int));
    out.write((char*)&(opts.nfeat),sizeof(int));
    out.write((char*)&(opts.ktype),sizeof(int));
    out.write((char*)&(opts.kpow),sizeof(int));
    out.write((char*)&(opts.kc1),sizeof(double));
    out.write((char*)&(opts.kc2),sizeof(double));
    out.write((char*)&(bias),sizeof(double));
    int itmp = opts.scale_kernel ? 1:0;
    out.write((char*)&(itmp),sizeof(int));
    for ( k = 0; k < nvecs; ++k ) {
        if ( status[k] < 0 ) continue;
        dtmp = y[k] * scal[k] * alfa[k];
        out.write((char*)&dtmp,sizeof(double));
    }
    for ( k = 0; k < nvecs ; ++k) {
        if (status[k]<0) continue;
        v = opts.vecs + k * nfeat;
        out.write((char*)v,nfeat*sizeof(double));
    }
    out.close();
    fprintf(stderr,"wrote model file\n");
    std::ofstream fout(opts.out.c_str());

    for ( k = 0; k < nvecs; ++k ) {
        fx = ( y[k] - grad[k] - bias );
        if ( fx > 0 ) {
            if ( y[k] > 0. ) ++ntp;
            else ++nfp;
        } else {
            if ( y[k] > 0. ) ++ntn;
            else ++nfn;
        }
        dtmp = y[k];
        fout.write((char*)&dtmp,sizeof(double));
        fout.write((char*)&fx,sizeof(double));
    }
    fout.close();
    analyze ( ntp, nfp, ntn, nfn );
}
