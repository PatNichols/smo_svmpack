#include "svm_classify.h"
#include "svm_data.h"

void svm_classify(svm_options_t *options) {
    int i;
    double sum;
    int ntp=0;
    int ntn=0;
    int nfp=0;
    int nfn=0;
    const double *vecs_data;
    const double *y_data;
    const double *vi;
    int nvecs_data;
    int nfeat_data;
    FILE *out;
    char * data_file;
    svm_model * model;
    svm_data * data;
#ifdef _OPENMP
    int narr[128][4];
    int nth_i[1];
    char *suffix;
    char *fname;
    int nsum;
    int fp,tp,fn,tn;
    int j;
    int xsz,sz;
    int tid,nth;
#endif    


    fprintf(stderr,"model file      = %s\n",options->model);
    data = svm_data_init(options->data);
    nfeat_data = data->nfeat;
    nvecs_data = data->nvecs;
    y_data = data->y;
    vecs_data = data->vecs;
    fprintf(stderr,"Classifying\n");
    fprintf(stderr,"# of data  vecs = %ld\n",nvecs_data);
    fprintf(stderr,"# of data  feat = %ld\n",nfeat_data);
    fprintf(stderr,"data file       = %s\n",options->data);
    fprintf(stderr,"out  file       = %s\n",options->out);
#ifdef _OPENMP
#pragma omp parallel private(fname,suffix,out,model,nth,tid,xsz,sz,j,tp,fp,tn,fn,i,vi,sum) shared(narr,nth_i)
  {
    fname = (char*)Malloc(32);
    suffix = (char*)Malloc(33);
    model = svm_model_init(options->model);
    nth = omp_get_num_threads();
    tid = omp_get_thread_num();
    if (tid==0) nth_i[0]=nth;
    xsz = nvecs_data%nth;
    sz = nvecs_data/nth;
    if (tid < xsz) {
        ++sz;
        j = sz*tid;
    }else{
        j = sz*tid + xsz;
    }
    sprintf(suffix,".%d\n",tid);
    fname = strcpy(fname,options->out);
    fname = strcat(fname,suffix);
    out = Fopen(fname,"w");
    tp = 0;
    fp = 0;
    tn = 0;
    fn = 0;
    for (i=0; i<sz; ++i) {
        vi = vecs_data + j * nfeat_data;
        sum  = svm_model_gensum(model,vi);
        if (sum >= 0.0) {
            /* positive */
            if (y_data[i]>0.) {
                ++tp;
            } else {
                ++fp;
            }
        } else {
            /* negative */
            if (y_data[i]>0.) {
                ++tn;
            } else {
                ++fn;
            }
        }
        fwrite((void*)&y_data[i],sizeof(double),1,out);
        fwrite((void*)&sum,sizeof(double),1,out);
        ++j;
    }
    fclose(out);
    narr[tid][0]=tp;
    narr[tid][1]=fp;
    narr[tid][2]=tn;
    narr[tid][3]=fn;
    svm_model_free(model);
    free(suffix);
    free(fname);
  }
  ntp = nfp = nfn = ntn = 0;
  nsum = 0;
  nth = nth_i[0];
  for (i=0;i<nth;++i) {
      nsum += narr[i][0] + narr[i][1] + narr[i][2] + narr[i][3];
      ntp += narr[i][0];
      nfp += narr[i][1];
      ntn += narr[i][2];
      nfn += narr[i][3];
      fprintf(stderr,"%3d %10d %10d %10d %10d\n",
          i,narr[i][0],narr[i][1],narr[i][2],narr[i][3]);
  }  
  fprintf(stderr,"nth = %d nsum = %d\n",nth,nsum);
#else
    model = svm_model_init(options->model);
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
    svm_model_free(model);
#endif
    analyze(ntp,ntn,nfp,nfn);
    svm_data_free(data);
}