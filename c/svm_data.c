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
        Fread((void*)&nvecs,sizeof(int),1,fp);
        Fread((void*)&nfeat,sizeof(int),1,fp);
        data->y = (double*)Malloc(nvecs*sizeof(double));
        data->vecs = (double*)Malloc(nvecs*nfeat*sizeof(double));
        Fread((void*)data->y,sizeof(double),nvecs,fp);
        Fread((void*)data->vecs,sizeof(double),nvecs*nfeat,fp);
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
    fprintf(stderr,"DATA \n");
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

