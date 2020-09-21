

#include "svm_classify.h"
#include "svm_data.h"

svm_data * svm_data_init(const char *data_file) {
    int nvecs,nfeat;
    int i,nt,nf;
    FILE *fp;
    svm_data * data = MALLOC_PTR(svm_data);
    
    if (strstr(data_file,".tdo")==0x0) {
        svm_data_read_libsvm(data,data_file);
    } else {
        fp = fopen(data_file,"r");
        fread((void*)&nvecs,1,sizeof(int),fp);
        fread((void*)&nfeat,1,sizeof(int),fp);
        data-> y = Malloc(nvecs*sizeof(double));
        data-> vecs = Malloc(nvecs*nfeat*sizeof(double));
        fread((void*)data->y,nvecs,sizeof(double),fp);
        fread((void*)data->vecs,nvecs*nfeat,sizeof(double),fp);
        data->nvecs = nvecs;
        data->nfeat = nfeat;
    }
    nt = 0;
    nf = 0;
    for (i=0; i<nvecs; ++i) {
        if (data->y[i] > 0.0) ++nt;
        else ++nf;
    }
    fprintf(stderr," data file # true = %ld # false = %ld\n",nt,nf);
    return data;
}

void svm_data_free(svm_data * data) {
    free(data->vecs);
    free(data->y);
    free(data);
}

void svm_data_read_libsvm(svm_data * data,const char *data_file) {
    int nvecs,nfeat;
    size_t buffer_size = MAX_LINE_SIZE;
    char *buffer = (char*)Malloc(MAX_LINE_SIZE);
    char **tokens = tokens_init();
    char *end;
    const char *delims = " :\n";
    size_t ntokens;
    size_t sz;
    int ivec,j;
    int index;
    double value;
    FILE *fp;

    nvecs = 0;
    nfeat = 0;
    fprintf(stderr,"Reading Data file %s\n",data_file);
    fp = Fopen(data_file,"r");
    while (!feof(fp)) {
        if (!getline(&buffer,&buffer_size,fp)) break;
        if (buffer[0]=='#') continue;
        ntokens = explode_string(buffer,delims,tokens);
        if (ntokens==0) break;
        if (ntokens%2) {
            index = ntokens-2;
            sz = atoi(tokens[index]);
            if (sz > nfeat) nfeat=(int)sz;
        } else {
            fprintf(stderr,"Format Error reading data file %s\n",data_file);
            exit(EXIT_FAILURE);
        }
        ++nvecs;
    }
    clearerr(fp);
    rewind(fp);
    fprintf(stderr,"# vecs = %ld # feat = %ld \n",nvecs,nfeat);
    data->nvecs =nvecs;
    data->nfeat =nfeat;
    data->y = (double*)Malloc(sizeof(double)*nvecs);
    data->vecs = (double*)Malloc(sizeof(double)*nvecs*nfeat);
    memset(data->vecs,0x0,sizeof(double)*nvecs*nfeat);
    ivec = 0;
    while (!feof(fp)) {
        if (!getline(&buffer,&buffer_size,fp)) break;
        if (buffer[0]=='#') continue;
        ntokens = explode_string(buffer,delims,tokens);
        if (ntokens==0) break;
        (data->y)[ivec] = strtod(tokens[0],&end);
        for (j=1; j<ntokens; j+=2) {
            index = atoi(tokens[j]);
            value = strtod(tokens[j+1],&end);
            (data->vecs)[ivec*nfeat + index-1] = value;
        }
        ++ivec;
    }
    fclose(fp);
    tokens_free(tokens);
    free(buffer);
}

void svm_model_free(svm_model *model) {
    free(model->vecs);
    free(model->yalfa);
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
    fread((void*)&nvecs,1,sizeof(int),fp);
    fread((void*)&nfeat,1,sizeof(int),fp);
    fread((void*)&(model->ktype),1,sizeof(int),fp);
    fread((void*)&(model->kpow),1,sizeof(int),fp);
    fread((void*)&(model->c1),1,sizeof(double),fp);
    fread((void*)&(model->c2),1,sizeof(double),fp);
    fread((void*)&(model->bias),1,sizeof(double),fp);
    fread((void*)&(model->kscal),1,sizeof(int),fp);
    vecs = (double*)Malloc(sizeof(double)*nvecs*nfeat);
    y = (double*)Malloc(sizeof(double)*nvecs);
    fread((void*)y,nvecs,sizeof(double),fp);
    fread((void*)vecs,nvecs*nfeat,sizeof(double),fp);
    fclose(fp);
    model->nfeat = nfeat;
    model->nvecs = nvecs;
    model->vecs = vecs;
    model->yalfa = y;
    if (model->ktype==2) model->c1 *= -1.0;
    nt = 0;
    nf = 0;
    for (i=0; i<nvecs; ++i) {
        if (y[i] > 0.0) ++nt;
        else ++nf;
    }
    fprintf(stderr," model file # true = %ld # false = %ld\n",nt,nf);
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
        for (i=0; i<nfeat; ++i) {
            sx += vp[i]*vp[i];
        }
        switch (model->ktype) {
        case 0:
            sx = (sx<tau) ? (1.0/sqrt(sx)):1.0;
            break;
        case 1:
            sx = model->c1 * sx + model->c2;
            sx = pow(sx,model->kpow);
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
        for (i=0; i<nvecs; ++i) {
            vi = model->vecs + i * nfeat;
            vsum = 0.0;
            for (j=0; j<nfeat; ++j) vsum += vp[j]*vi[j];
            sum += vsum * (model->yalfa)[i];
        }
        return sum*sx-model->bias;
    case 1:
        for (i=0; i<nvecs; ++i) {
            vi = model->vecs + i * nfeat;
            vsum = 0.0;
            for (j=0; j<nfeat; ++j) vsum += vp[j]*vi[j];
            vsum = pow(model->c1 * vsum + model->c2,model->kpow);
            sum += vsum * (model->yalfa)[i];
        }
        return sum*sx-model->bias;
    case 2:
        sum = 0.0;
        for (i=0; i<nvecs; ++i) {
            vi = model->vecs + i * nfeat;
            vsum = 0.0;
            for (j=0; j<nfeat; ++j) {
                t = vp[j]-vi[j];
                vsum += t*t;
            }
            vsum = exp(model->c1*vsum);
            sum += vsum * (model->yalfa)[i];
        }
        return sum-(model->bias);
    case 3:
        for (i=0; i<nvecs; ++i) {
            vi = model->vecs + i * nfeat;
            vsum = 0.0;
            for (j=0; j<nfeat; ++j) vsum += vp[j]*vi[j];
            vsum = tanh(model->c1 * vsum + model->c2);
            sum += vsum * (model->yalfa)[i];
        }
        return sum*sx-model->bias;
    }
}

void svm_classify(svm_options_t *options) {
    int i;
    double sum;
    int ntp=0;
    int ntn=0;
    int nfp=0;
    int nfn=0;
    double *vecs_data;
    double *y_data;
    const double *vi;
    int nvecs_data;
    int nfeat_data;
    FILE *fp;
    FILE *out;
    char * data_file;
    svm_model * model;
    svm_data * data;

    fprintf(stderr,"model file      = %s\n",options->model);
    model = svm_model_init(options->model);
    data = svm_data_init(options->data);
    nfeat_data = data->nfeat;
    nvecs_data = data->nvecs;
    y_data = data->y;
    vecs_data = data->vecs;
    fprintf(stderr,"Classifying\n");
    fprintf(stderr,"# of model vecs = %ld\n",model->nvecs);
    fprintf(stderr,"# of model feat = %ld\n",model->nfeat);
    fprintf(stderr,"# of data  vecs = %ld\n",nvecs_data);
    fprintf(stderr,"# of data  feat = %ld\n",nfeat_data);
    fprintf(stderr,"data file       = %s\n",options->data);
    fprintf(stderr,"out  file       = %s\n",options->out);

#ifdef _OPENMP

#pragma omp parallel for private(i,vi,sum) reduction(+:ntp) reduction(+:nfp) reduction(+:ntn) reduction(+:nfn)
    for (i=0; i<nvecs_data; ++i) {
        vi = vecs_data + i * nfeat_data;
        sum  = svm_model_gensum(model,vi);
        if (sum >= 0.0) {
            /* positive */
            if (y_data[i]>0.) {
                ++ntp;
            } else {
                ++nfp;
            }
        } else {
            /* negative */
            if (y_data[i]>0.) {
                ++ntn;
            } else {
                ++nfn;
            }
        }
    }
#else
    ntp = 0;
    nfp = 0;
    ntn = 0;
    nfn = 0;
    out = Fopen(options->out,"w");
    for (i=0; i<nvecs_data; ++i) {
        vi = vecs_data + i * nfeat_data;
        sum  = svm_model_gensum(model,vi);
        if (sum >= 0.0) {
            /* positive */
            if (y_data[i]>=0.) {
                ++ntp;
            } else {
                ++nfp;
            }
        } else {
            /* negative */
            if (y_data[i]>=0.) {
                ++ntn;
            } else {
                ++nfn;
            }
        }
        fwrite((void*)&y_data[i],1,sizeof(double),out);
        fwrite((void*)&sum,1,sizeof(double),out);
    }
    fclose(out);
#endif
    analyze(ntp,ntn,nfp,nfn);
    svm_data_free(data);
    svm_model_free(model);
}