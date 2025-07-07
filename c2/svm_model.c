#include "svm_model.h"
#include <math.h>
#include "stopwatch.h"

static inline double model_powi(double x,int m) {
    double w,y,z;
    if (m<0) {
        m = -m;
        x = 1./x;
    }
    switch(m) {
    case 0:
        return 1.0;
    case 1:
        return x;
    case 2:
        return x*x;
    case 3:
        return x*x*x;
    case 4:
        y =x*x;
        return y*y;
    case 5:
        y = x*x;
        return y*y*x;
    case 6:
        y = x*x;
        return y*y*y;
    case 7:
        y = x*x;
        return y*y*y*x;
    case 8:
        y= x*x;
        y= y*y;
        return y*y;
    case 9:
        y = x*x;
        y = y*y;
        return y*y*x;
    case 10:
        y = x*x;
        z = y*y;
        z = y*y;
        return z*y;
    case 11:
        y = x*x;
        z = y*y;
        z = y*y;
        return z*y*x;
    case 12:
        y = x*x;
        z = y*y*y;
        return z*z;
    default:
        y = x*x;
        z = y*y*y;
        z = z*z;
        return z*model_powi(x,(m-12));
    }
    return 1.e300;
}

void svm_model_free(svm_model *model) {
    free(model->vecs);
    free(model->yalf);
    free(model);
}

svm_model * svm_model_init(const char * model_file,int tid) {
    double *vecs;
    double *yalf;
    int nvecs;
    int nfeat;
    int i,nt,nf,itmp;
    size_t vsize;
    svm_model* model = (svm_model*)Malloc(sizeof(svm_model));
    FILE *in=Fopen(model_file,"r");
    Fread((void *)&nvecs,sizeof(int),1,in);
    Fread((void *)&nfeat,sizeof(int),1,in);
    Fread((void *)&(model->ktype),sizeof(int),1,in);
    Fread((void *)&(model->kpow),sizeof(int),1,in);
    Fread((void *)&(model->c1),sizeof(double),1,in);
    Fread((void *)&(model->c2),sizeof(double),1,in);
    Fread((void *)&(model->bias),sizeof(double),1,in);
    Fread((void *)&itmp,sizeof(int),1,in);
    if (itmp==0) model->kscal = 0;
    else model->kscal = 1;
    vsize = nvecs * nfeat;
    yalf = (double*)Malloc(sizeof(double)*nvecs);
    vecs = (double*)Malloc(sizeof(double)*vsize);
    Fread((void*)yalf,sizeof(double),nvecs,in);
    Fread((void*)vecs,sizeof(double),vsize,in);
    fclose(in);
    model->nvecs = nvecs;
    model->nfeat = nfeat;
    model->yalf = yalf;
    model->vecs = vecs;
    if (model->ktype<0 || model->ktype>3) {
        fprintf(stderr,"unknown kernel type in model init %d\n",model->ktype);
        exit(EXIT_FAILURE);
    }
    if (model->ktype == 2) model->c1 *= -1.;

    if (tid==0) {
        nt = 0;
        nf = 0;
        for (i=0; i<nvecs; ++i) {
            if ( (model->yalf)[i]>=0.) ++nt;
            else ++nf;
        }
        fprintf(stderr,"svm model\n");
        fprintf(stderr,"model file     = %s\n",model_file);
        fprintf(stderr,"nvecs          = %d\n",model->nvecs);
        fprintf(stderr,"nfeat          = %d\n",model->nfeat);
        fprintf(stderr,"kernel type    = %d\n",model->ktype);
        fprintf(stderr,"kernel power   = %d\n",model->nvecs);
        fprintf(stderr,"c1             = %lg\n",model->c1);
        fprintf(stderr,"c2             = %lg\n",model->c2);
        fprintf(stderr,"kernel scale   = %d\n",model->kscal);
        fprintf(stderr,"bias           = %lg\n",model->bias);
        fprintf(stderr,"n positive     = %d\n",nt);
        fprintf(stderr,"n negative     = %d\n",nf);
    }
    return model;
}

double svm_model_scale_factor(svm_model* model,const double *vp) {
    double sum;
    double t;
    const double eps=1.e-14;

    double scal = 1.0;
    if (model->ktype==2 || model->kscal==0) return 1.0;
    switch (model->ktype) {
    case 0:
        sum = svm_dot(vp,vp,model->nfeat);
        break;
    case 1:
        sum = svm_dot(vp,vp,model->nfeat)*(model->c1) + model->c2;
        sum = model_powi(sum,model->kpow);
        break;
    case 2:
        sum = 1.0;
        break;
    case 3:
        sum = svm_dot(vp,vp,model->nfeat)*(model->c1) + model->c2;
        sum = tanh(sum);
        break;
    default:
        break;
    }
    sum = sqrt(sum);
    if (sum > eps) return 1./sum;
    else return 1.0;
}

double svm_model_kernel_sum(svm_model* model,const double *vp) {
    int i;
    double scal_vp;
    double sum;
    int nvecs = model->nvecs;
    int nfeat = model->nfeat;
    const double *vecs = model->vecs;
    const double *yalf = model->yalf;
    const double bias = model->bias;
    const double c1 = model->c1;
    const double c2 = model->c2;
    int kpow = model->kpow;
    switch(model->ktype) {
    case 0:
        scal_vp = svm_model_scale_factor(model,vp);
        sum = 0.0;
        for (i=0; i<nvecs; ++i) {
            sum += svm_dot(vecs+i*nfeat,vp,nfeat)*yalf[i];
        }
        sum = sum*scal_vp - bias;
        return sum;
    case 1:
        scal_vp = svm_model_scale_factor(model,vp);
        sum = 0.0;
        for (i=0; i<nvecs; ++i) {
            sum += model_powi(c1*svm_dot(vecs+i*nfeat,vp,nfeat)+c2,kpow)*yalf[i];
        }
        sum = sum*scal_vp - bias;
        return sum;
    case 2:
        sum = 0.0;
        for (i=0; i<nvecs; ++i) {
            sum += exp(c1*svm_diff_nrm2(vecs+i*nfeat,vp,nfeat))*yalf[i];
        }
        sum -= bias;
        return sum;
    case 3:
        scal_vp = svm_model_scale_factor(model,vp);
        sum = 0.0;
        for (i=0; i<nvecs; ++i) {
            sum += tanh(c1*svm_dot(vecs+i*nfeat,vp,nfeat)+c2)*yalf[i];
        }
        sum = sum*scal_vp - bias;
        return sum;
    }
}

void svm_model_classify(svm_options_t *options) {
    svm_model * model;
    svm_data * data;
    double fx;
    int ntp,nfp,ntn,nfn;
    int i;
    FILE *fout;
    double *yx;
    double *vecs;
    int nfeat;
    int nvecs;
    stopwatch_t rd_timer;
    stopwatch_t t_timer;
#ifdef _OPENMP
    int nth,tid,sz,xsz,k,ntot;
    int nth_i[1];
    char *file_name;
    char suffix[15];
    int narr[128][4];
#endif
    stopwatch_clear(&rd_timer);
    stopwatch_clear(&t_timer);
    stopwatch_start(&rd_timer);
    data = svm_data_init(options->data);
    stopwatch_stop(&rd_timer);
    fprintf(stderr,"time to read data = %lg seconds\n",stopwatch_time(&rd_timer));
    stopwatch_start(&t_timer);
    yx = data->y;
    vecs = data->vecs;
    nvecs = data->nvecs;
    nfeat = data->nfeat;
    ntp = 0;
    nfp = 0;
    nfn = 0;
    ntn = 0;
#ifdef _OPENMP
    #pragma omp parallel private(nth,tid,sz,k,i,file_name,suffix,fout,fx,model,ntp,nfp,ntn,nfn) shared(nth_i)
    {
        nth = omp_get_num_threads();
        if (tid==0) *nth_i = nth;
        tid = omp_get_thread_num();
        ntp = 0;
        ntn = 0;
        nfp = 0;
        nfn = 0;
        narr[tid][0]= 0;
        narr[tid][1]= 0;
        narr[tid][2]= 0;
        narr[tid][3]= 0;
        model = svm_model_init(options->model,tid);
        sz = nvecs/nth;
        xsz = nvecs%nth;
        if (tid < xsz) {
            sz +=1;
            k = sz*tid;
        } else {
            k = sz*tid + xsz;
        }
        file_name = (char*)Malloc(128);
        file_name = strncpy(file_name,options->out,115);
        snprintf(suffix,11,"%d",tid);
        strncat(file_name,suffix,11);
        fout = Fopen(options->out,"w");
        for (i=0; i<sz; ++i,++k) {
            fx = svm_model_kernel_sum(model,vecs+k*nfeat);
            if (fx >= 0.0) {
                if (yx[k]>=0.0) ++ntp;
                else ++nfp;
            } else {
                if (yx[k]>=0.0) ++nfn;
                else ++ntn;
            }
            fwrite((void*)&yx[k],sizeof(double),1,fout);
            fwrite((void*)&fx,sizeof(double),1,fout);
        }
        fclose(fout);
        free(file_name);
        svm_model_free(model);
        narr[tid][0]= ntp;
        narr[tid][1]= ntn;
        narr[tid][2]= nfn;
        narr[tid][3]= nfp;
    }
    ntp = 0;
    ntn = 0;
    nfp = 0;
    ntn = 0;
    ntot = 0;
    nth = *nth_i;
    for (i=0; i<nth; ++i) {
        ntp += narr[i][0];
        ntn += narr[i][1];
        nfn += narr[i][2];
        nfp += narr[i][3];
        ntot += narr[i][0] + narr[i][1] + narr[i][2] + narr[i][3];
    }
    stopwatch_stop(&t_timer);
    fprintf(stderr,"time for classification = %lg seconds\n",stopwatch_time(&t_timer));
#else
    stopwatch_start(&t_timer);
    model = svm_model_init(options->model,0);
    fout = Fopen(options->out,"w");
    for (i=0; i<nvecs; ++i) {
        fx = svm_model_kernel_sum(model,vecs+i*nfeat);
        if (fx >= 0.0) {
            if (yx[i]>=0.0) ++ntp;
            else ++nfp;
        } else {
            if (yx[i]>=0.0) ++nfn;
            else ++ntn;
        }
        fwrite((void*)&yx[i],sizeof(double),1,fout);
        fwrite((void*)&fx,sizeof(double),1,fout);
    }
    fclose(fout);
    stopwatch_stop(&t_timer);
    fprintf(stderr,"time for classification = %lg seconds\n",stopwatch_time(&t_timer));
    svm_model_free(model);
    ntot = nfp + nfn + ntp + ntn;
#endif
    fprintf(stderr,"nt = %d\n",ntot);
    fprintf(stderr,"ntp = %d\n",ntp);
    fprintf(stderr,"ntn = %d\n",ntn);
    fprintf(stderr,"nfp = %d\n",nfp);
    fprintf(stderr,"nfn = %d\n",nfn);
    analyze(ntp,ntn,nfp,nfn);
    svm_data_free(data);
}
