#include "svm_data.h"

svm_data * svm_data_init(const char *data_file) {
    int nvecs,nfeat,nt,nf;
    int i;
    FILE *fp;
    svm_data * data = MALLOC_PTR(svm_data);
    
    if (strstr(data_file,".tdo")==0x0) {
        read_libsvm_file(data_file,&nvecs,&nfeat,&(data->y),&(data->vecs));
        data->nvecs = nvecs;
        data->nfeat = nfeat;
    }else{
        fp = Fopen(data_file,"r");
        fread((void*)&nvecs,sizeof(int),1,fp);
        fread((void*)&nfeat,sizeof(int),1,fp);
        data-> y = (double*)Malloc(nvecs*sizeof(double));
        data-> vecs = (double*)Malloc(nvecs*nfeat*sizeof(double));
        fread((void*)data->y,sizeof(double),nvecs,fp);
        fread((void*)data->vecs,sizeof(double),nvecs*nfeat,fp);
        fclose(fp);
        data->nvecs = nvecs;
        data->nfeat = nfeat;
    }
    nt = 0;
    nf = 0;
    for (i=0;i<nvecs;++i) {
        if (data->y[i] > 0.0) ++nt;
        else ++nf;
    }

    fprintf(stderr,"data file     = %s\n",data_file);
    fprintf(stderr,"# of vecs     = %d\n",nvecs);
    fprintf(stderr,"# of features = %d\n",nfeat);
    fprintf(stderr,"# true        = %d\n",nt);
    fprintf(stderr,"# false       = %d\n\n",nf);
    return data;
}

void svm_data_free(svm_data* data) {
    free(data->vecs);
    free(data->y);
    free(data);
}

void svm_model_free(svm_model *model) {
    free(model->vecs);
    free(model->y);
    free(model);
}


svm_model * svm_model_init(const char *model_file) {
    double *y;
    double *vecs;
    int nvecs;
    int nfeat;
    int nf,nt,i;
    size_t nrd;
    /// classifying task read in model file
    svm_model *model = MALLOC_PTR(svm_model);
    FILE *fp = Fopen(model_file,"r");
    fread((void*)&nvecs,sizeof(int),1,fp);
    fread((void*)&nfeat,sizeof(int),1,fp);
    fread((void*)&(model->ktype),sizeof(int),1,fp);
    fread((void*)&(model->kpow),sizeof(int),1,fp);
    fread((void*)&(model->c1),sizeof(double),1,fp);
    fread((void*)&(model->c2),sizeof(double),1,fp);
    fread((void*)&(model->bias),sizeof(double),1,fp);
    fread((void*)&(model->kscal),sizeof(int),1,fp);
    vecs = (double*)Malloc(sizeof(double)*nvecs*nfeat);
    y = (double*)Malloc(sizeof(double)*nvecs);
    fread((void*)y,sizeof(double),nvecs,fp);
    fread((void*)vecs,sizeof(double),nvecs*nfeat,fp);
    fclose(fp);
    model->nfeat = nfeat;
    model->nvecs = nvecs;
    model->vecs = vecs;
    model->y = y;
    if (model->ktype==2) model->c1 = -(model->c1);
    nt = 0;
    nf = 0;
    for (i=0;i<nvecs;++i) {
        if (y[i] > 0.0) ++nt;
        else ++nf;
    }

    fprintf(stderr,"model file = %s\n",model_file);
    fprintf(stderr,"# of vecs     = %d\n",nvecs);
    fprintf(stderr,"# of features = %d\n",nfeat);
    fprintf(stderr,"# true        = %d\n",nt);
    fprintf(stderr,"# false       = %d\n",nf);
    fprintf(stderr,"kernel type   = %d\n",model->ktype);
    fprintf(stderr,"kernel power  = %d\n",model->kpow);
    fprintf(stderr,"kernel scale  = %lg\n",model->kscal);
    fprintf(stderr,"kernel c1     = %lg\n",model->c1);
    fprintf(stderr,"kernel c2     = %lg\n",model->c2);
    fprintf(stderr,"bias          = %lg\n\n",model->bias);
    return model;
}

double svm_model_gensum(svm_model *model,const double *vp)
{    
    int i,j;
    int nvecs = model->nvecs;
    int nfeat = model->nfeat;
    double sx = 1.;
    double vsum,sum,t;
    const double tau=1.e-14;
    const double *vi;
    if (model->kscal && model->ktype!=2) {
        sx = 0.0;
        for (i=0;i<nfeat;++i) {    
            sx += vp[i]*vp[i];
        }
        switch (model->ktype) {
        case 0:
            sx = (sx<tau) ? (1.0/sqrt(sx)):1.0;
            break;
        case 1:
            sx = model->c1 * sx + model->c2;
            sx = dpowi(sx,model->kpow);
            sx = (sx<tau) ? (1.0/sqrt(sx)):1.0;
            break;
        case 2:
            sx = 1.0;
            break;
        case 3:
            sx = model->c1 * sx + model->c2;
            sx = tanh(sx);
            sx = (sx<tau) ? (1.0/sqrt(sx)):1.0;
            break;
        }
    }
    switch(model->ktype) {
    case 0:
        for (i=0;i<nvecs;++i) {
            vi = model->vecs + i * nfeat;
            vsum = 0.0;
            for (j=0;j<nfeat;++j) vsum += vp[j]*vi[j];
            sum += vsum * (model->y)[i];
        }
        return sum*sx-model->bias;
    case 1:
        for (i=0;i<nvecs;++i) {
            vi = model->vecs + i * nfeat;
            vsum = 0.0;
            for (j=0;j<nfeat;++j) vsum += vp[j]*vi[j];
            vsum = dpowi(model->c1 * vsum + model->c2,model->kpow);
            sum += vsum * (model->y)[i];
        }
        return sum*sx-model->bias;
    case 2:
        for (i=0;i<nvecs;++i) {
            vi = model->vecs + i * nfeat;
            vsum = 0.0;
            for (j=0;j<nfeat;++j) {
                t = vp[j]-vi[j];
                vsum += t*t;
            }
            vsum = exp(model->c1*vsum);
            sum += vsum * (model->y)[i];
        }
        return sum-model->bias;
    case 3:
        for (i=0;i<nvecs;++i) {
            vi = model->vecs + i * nfeat;
            vsum = 0.0;
            for (j=0;j<nfeat;++j) vsum += vp[j]*vi[j];
            vsum = tanh(model->c1 * vsum + model->c2);
            sum += vsum * (model->y)[i];
        }
        return sum*sx-model->bias;
    }
}

