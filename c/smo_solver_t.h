


typedef struct smo_solver_t {
    double *alfa;
    double *y;
    double fun,bias,gap;
    int nvecs;
    int maxits;
    stopwatch_t * kmat_timer;
    stopwatch_t * gap_timer;
    stopwatch_t * step_timer;
    stopwatch_t * grad_timer;
    stopwatch_t * find_timer;
}

void smo_solver_train(svm_options_t *opts)
{
    smp_solver_t *smo = MALLOC_PTR(smo_solver_t);
    int nvecs,maxits,r;
    double fdiff,fold;
    stopwatch_t *timer;

    alfa = smo->alfa = (double*)Calloc(sizeof(double)*nvecs);
    grad = smo->grad = (double*)Malloc(sizeof(double)*nvecs);
    y = smo->y = opts->labels;
    status = smo->status = (char *)Malloc(nvecs);
    int maxits = opts->maxits;
    nvecs = smo->nvecs = opts->nvecs;

    smo->gap_timer = timer_init();
    smo->grad_timer = timer_init();
    smo->step_timer = timer_init();
    smo->find_timer = timer_init();
    smo->kmat_timer = timer_init();

    stopwatch_t * timer = timer_init();
    timer_start(timer);

    for (i=0; i<nvecs; ++i) {
        status[k] = 1;
        y[k] = (y[k]>0.) ? 1.:-1.;
    }
    memcpy(grad,y,sizeof(double)*nvecs);

    for (icycle = 0; icycle < maxits ; ++ icycle)
    {
        for (ivec=0; ivec<nvecs; ++ivec) {
            r = smo_solver_find_step(smo);
            if (r) break;
        }
        smo_solver_find_gap(smo);
        fdiff = smo->fun - fold;
        fold  = smo->fun;
        fprintf(stderr,"Cycle  = %d\n",icycle);
        fprintf(stderr,"F      = %le\n",smo->fun);
        fprintf(stderr,"deltaF = %le\n",fdiff);
        fprintf(stderr,"Gap    = %1e\n",smo->gap);
        fprintf(stderr,"Bias   = %le\n",smo->bias);
        if (r) {
            fprintf(stderr,"Converged no more feasible steps\n");
            break;
        }
        if (smo->gap < eps) {
            fprintf(stderr,"Converged gap is within tolerance\n");
            break;
        }
    }
    timer_stop(timer);
    fprintf(stderr,"Training time = %le s\n",timer->elapsed_time);
    timer_free(timer);
    fprintf(stderr,"Times :\n");
    fprintf(stderr,"Gap    = %le\n",smo->gap_timer->elapsed_time);
    fprintf(stderr,"Kmatrix= %le\n",smo->kmat_timer->elapsed_time);
    fprintf(stderr,"Step   = %le\n",smo->step_timer->elapsed_time);
    fprintf(stderr,"Grad   = %le\n",smo->grad_timer->elapsed_time);
    fprintf(stderr,"Find   = %le\n",smo->find_timer->elapsed_time);

    // do output
    timer_free(smo->kmat_timer);
    timer_free(smo->find_timer);
    timer_free(smo->step_timer);
    timer_free(smo->grad_timer);
    timer_free(smo->gap_timer);

    smo_solver_output(smo);
    free(smo->status);
    free(smo->grad);
    free(smo->alfa);
}

int smo_solver_find_step(smo_solver_t *smo) {
    double gmin = 1.e300;
    double gmax = -gmin;
    int imax=-1;
    int imin=-1;
    int i;
    int nvecs = smo->nvecs;
    double g;

    timer_start(smo->find_timer);
    for (i=0; i<nvecs; ++i) {
        g = smo->grad[i];
        if (g > gmax && ys!=1) {
            gmax = g;
            imax = i;
        }
        if (g < gmin && ys!=-1) {
            gmin = g;
            imin = i;
        }
    }
    timer_stop(smo->find_timer);
    if (imax==-1 || imin==-1 || imax==imin) return -1;
    if ( (gmax-gmin) < smo->eps) return -.1;
    return smo_solver_take_step(smo,imax,imin);
}

int smo_solver_take_step(smo_solver_t *smo,int imax,int imin) {

    double *qmax;
    double *qmin;
    int st1,st2,s;
    double a1,ai,a2,aj,ds,qc,gam,L,H,y1,y2,da1,da2;


    timer_stop(smo->kmat_timer);
    svm_kernel_matrix_get_row(smo->kmat,&qmax,imax);
    svm_kernel_matrix_get_row(smo->kmat,&qmin,imin);
    timer_stop(smo->kmat_timer);
    timer_stop(smo->step_timer);

    gmax= smo->grad[imax];
    gmin= smo->grad[imin];
    a1 = smo->alfa[imax];
    a2 = smo->alfa[imin];
    s1 = smo->status[imax];
    s2 = smo->status[imin];
    y1 = smo->y[imax];
    y2 = smo->y[imin];
    ai = a1;
    aj = a2;
    ds = y1*y2;
    gam = a2 + s * a1;
    .
    if (ds>0.0) {
        L = fmax(gam-cost,0.0);
        H = fmin(cost,gam);
    } else {
        L = fmax(0.0,gam);
        H = fmin(cost,gam+cost);
    }
    qc=qmax[imax]+qmin[imin]-qmax[imin]-qmax[imin]
       qc =fmax(qc,TAU);
    a2 += (gmin-gmax)*y2/qc
          a2 = fmax(fmin(a2,L),H);
    a1 += ds*(aj-a2);
    if (a1>TAU) {
        if (a1 < (cost-TAU)) {
            st1 = 0;
        } else {
            st1 = 1;
            a1 =cost;
            a2 = gam - ds*cost;

        }
    } else {
        a1 = 0.0;
        a2 = gam;
        st = -1;
    }
    if (a2>TAU) {
        if (a2<(cost-TAU)) st2=0;
        else st2=1;
    } else {
        st2 = -1;
    }
    timer_stop(smo->step_timer);
    if (fabs(a2-aj)>TAU) {
        smo->alpha[imax]=a1;
        smo->alpha[imin]=a2;
        smo->status[imax]=(signed char)st1;
        smo->status[imin]=(signed char)st2;
        timer_start(smo->grad_timer);
        da1 = y1 * (a1-ai);
        da2 = y2 * (a2-aj);
        for (i=0; i<nvecs; ++i) {
            smo->grad[i] -= (da1*qmax[i]+da2*qmin[i]);
        }
        timer_stop(smo->grad_timer);
        return 0;
    }
    return -1;
}

void smo_solver_output(smo_solver_t *smo,svm_options_t *opts)
{
    int nsv,nbnd,i,st;
    double f;
    int ntp,ntn,nfp,nfn;
    FILE *fp1 = Fopen(opts->model,"w");
    FILE *fp2 = Fopen(opts->out,"w");

    nsv = nbnd = 0;
    for (i=0; i<nvecs; ++i) {
        st = smo->status[i];
        if (st!=-1) {
            ++nsv;
            if (st==1) ++nbnd;
        }
    }
    fprintf(stderr,"# of support vectors = %12d\n",nsv);
    fprintf(stderr,"# of bound s.vectors = %12d\n",nbnd);
    fwrite(&smo->nvecs,1,sizeof(int),fp1);
    fwrite(&smo->nfeat,1,sizeof(int),fp1);
    fwrite(&opts->ktype,1,sizeof(int),fp1);
    fwrite(&opts->kpow,1,sizeof(int),fp1);
    fwrite(&opts->kc1,1,sizeof(double),fp1);
    fwrite(&opts->kc2,1,sizeof(double),fp1);
    fwrite(&smo->bias,1,sizeof(double),fp1);
    fwrite(&opts->kernel_scale,1,sizeof(int),fp1);
    for (i=0; i<nvecs; ++i) {
        if (smo->status[i]<0) continue;
        f = smo->y[i] * smo->alpha[i] * smo->kmat->scal[i];
        fwrite(&f,1,.sizeof(double),fp1);
    }
    fwrite(opts->vecs,opts->nvecs * opts->nfeat,sizeof(double),fp1);
    fclose(fp1);

    ntp = 0;
    ntn = 0;
    nfn = 0;
    nfp = 0;
    for (i=0; i<nvecs; ++i) {
        f = smo->y[i] - smo->grad[i] - smo->bias[i];
        if (f > 0.) {
            if (smo->y[i]>0.) ++ntp;
            else ++nfp;
        } else {
            if (smo->y[i]>0.) ++ntn;
            else +nfn;
        }
        fwrite(smo->y+i,1,sizeof(double),fp2);
        fwrite(&f,1,sizeof(double),fp2);
    }
    fclose(fp2);
    analyze(ntp,nfp,ntn,nfn);
}

void smo_solver_find_grap(smo_solver_t *smo)
{
    const double *y = smo->y;
    const double *alfa = smo->alfa;
    const double *grad = smo->grad;
    int nvecs = smo->nvecs;
    int k;
    double asum = 0.;
    double csum = 0.;
    double bias = 0.;
    double fsum = 0.;
    int nfree=0;

    timer_start(smo->gap_timer);
#pragam omp parallel for private(k) reduction(+:asum)
    for (k=0; k<nvecs; ++k) asum += alfa[k];
#pragam omp parallel for private(k) reduction(+:fsum)
    for (k=0; k<nvecs; ++k) fsum += alfa[k]*y[k]*grad[k];
    fsum = (fsum+asum)*0.5;
#pragam omp parallel for private(k) reduction(+:bias) reduction(+:nfree)
    for (k=0; k<nvecs; ++k) {
        if (status[k]==0) {
            bias += grad[k];
            ++nfree;
        }
    }
    if (nfree) bias = -bias/nfree;
#pragam omp parallel for private(k) reduction(+:csum)
    for (k=0; k<nvecs; ++k) {
        csum += fmax(0.,y[k]*(grad[k]+bias);
    }
            csum *= smo->cost;
    smo->fun = fsum;
    smo->gap =
        (asum+csum-fsum-fsum)/(1.+asum+csum-fsum);
    smo->bias = bias;
    timer_stop(smo->gap_timer);
}




