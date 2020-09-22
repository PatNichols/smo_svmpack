#include "svm_options.h"



inline void svm_options_read_model_file(svm_options_t* opts)
{
    double *y;
    double *vecs;
    int nvecs;
    int nfeat;
    int nf,nt,i;
    size_t nrd;
    /// classifying task read in model file
    FILE *fp = Fopen(opts->model,"r");
    Fread((void*)&nvecs,sizeof(int),1,fp);
    Fread((void*)&nfeat,sizeof(int),1,fp);
    Fread((void*)&(opts->ktype),sizeof(int),1,fp);
    Fread((void*)&(opts->kpow),sizeof(int),1,fp);
    Fread((void*)&(opts->kc1),sizeof(double),1,fp);
    Fread((void*)&(opts->kc2),sizeof(double),1,fp);
    Fread((void*)&(opts->bias),sizeof(double),1,fp);
    Fread((void*)&(opts->scale_kernel),sizeof(int),1,fp);
    vecs = (double*)Malloc(sizeof(double)*nvecs*nfeat);
    y = (double*)Malloc(sizeof(double)*nvecs);
    Fread((void*)y,sizeof(double),nvecs,fp);
    Fread((void*)vecs,sizeof(double),nvecs*nfeat,fp);
    fclose(fp);
    opts->nfeat = nfeat;
    opts->nvecs = nvecs;
    opts->vecs = vecs;
    opts->y = y;
    nt = 0;
    nf = 0;
    for (i=0;i<nvecs;++i) {
        if (y[i] > 0.0) ++nt;
        else ++nf;
    }
    fprintf(stderr," model file # true = %ld # false = %ld\n",nt,nf);
    return;
}

inline void svm_options_translate(svm_options_t* options) {
    char *pstr = strstr(options->data,".tdo");

    fprintf(stderr,"translating\n");
    fprintf(stderr,"data file is %s\n",options->data);
    fprintf(stderr,"out file is %s\n",options->out);
    svm_options_read_data_file(options);
    if (strstr(options->data,".tdo")!=0x0) {
        /* data file is a tdo file, write out libsvm file */
        write_libsvm_file(options->out,options->nvecs,options->nfeat,options->y,options->vecs);
    }else{
        write_tdo_file(options->out,options->nvecs,options->nfeat,options->y,options->vecs);
    }
}

int find_substring(const char *str,const char *sub)
{
    char * s = strstr(str,sub);
    if (s) return (s-str);
    return -1;
}

inline void svm_options_read_data_file(svm_options_t* opts) {
    int nt;
    int nf;
    int i;
    char * px= strstr(opts->data,".tdo");

    fprintf(stderr,"reading data file %s\n",opts->data);
    if (px) {
        read_tdo_file(opts->data,&(opts->nvecs),&(opts->nfeat),&(opts->y),&(opts->vecs));
    }else{
        read_libsvm_file(opts->data,&(opts->nvecs),&(opts->nfeat),&(opts->y),&(opts->vecs));
    }
    fprintf(stderr,"read data file %s\n",opts->data);
    nt = 0;
    nf = 0;
    for (int i=0;i<opts->nvecs ;++i) {
        if (opts->y[i]>0.0) ++nt;
        else ++nf;
    }
    fprintf(stderr,"# true = %ld #false = %ld\n",nt,nf);
}

inline svm_options_t * svm_options_init(int argc,char **argv)
{
    int task;
    char *end;
    char *task_str;
    svm_options_t * opts;
    program_options_t * popts = program_options_init(argc,argv);

    program_options_insert(popts,"data","input data file","svm.in");
    program_options_insert(popts,"model","model file name","svm.model");
    program_options_insert(popts,"out","output file name","svm.out");
    program_options_insert(popts,"task","task to perform (train,classify,translate)","train");
    program_options_insert(popts,"kernel_type","kernel_function type:\n   0)dot product\n   1)polynomial dot product\
        \n   2)radial basis function\n   3)logistic function","2");
    program_options_insert(popts,"kernel_power","kernel power for type 1 kernel function","2");
    program_options_insert(popts,"kernel_cof1","first parameter for kernel function","-1.0");
    program_options_insert(popts,"kernel_cof2","second parameter for kernel function","0.0");
    program_options_insert(popts,"cost","cost parameter for soft margin training","1.0");
    program_options_insert(popts,"eps","convergence parameter for training","1.e-12");
    program_options_insert(popts,"nthreads","number of threads","0");
    program_options_insert(popts,"cache_size","number of rows to cache","-1");
    program_options_insert(popts,"scale","scale kernel so diagonal elements are 1","true");
    program_options_insert(popts,"config","config file for options",0x0);
    program_options_insert(popts,"max_iterations","max # of training cycles","0");
    program_options_parse_command_line(popts,argc,argv);
    if (program_options_has_value(popts,"config"))
        program_options_parse_config_file(popts,program_options_get_value(popts,"config"));

    opts = MALLOC_PTR(svm_options_t);
    opts->data = strdup(program_options_get_value(popts,"data"));
    opts->model = strdup(program_options_get_value(popts,"model"));
    opts->out = strdup(program_options_get_value(popts,"out"));
    task_str = strdup(program_options_get_value(popts,"task"));
    opts->ktype = atoi(program_options_get_value(popts,"kernel_type"));
    opts->kpow = atoi(program_options_get_value(popts,"kernel_power"));
    opts->kc1 = strtod(program_options_get_value(popts,"kernel_cof1"),&end);
    opts->kc2 = strtod(program_options_get_value(popts,"kernel_cof2"),&end);
    opts->cost = strtod(program_options_get_value(popts,"cost"),&end);
    opts->eps = strtod(program_options_get_value(popts,"eps"),&end);
    opts->nths = atoi(program_options_get_value(popts,"nthreads"));
    opts->csize = atoi(program_options_get_value(popts,"cache_size"));
    opts->scale_kernel = parse_bool(program_options_get_value(popts,"scale"));
    opts->max_its = atoi(program_options_get_value(popts,"max_iterations"));

    if (strcmp(task_str,"train")==0) {
        opts->task=0;
        svm_options_read_data_file(opts);
        if (opts->csize == -1) opts->csize = opts->nvecs/6;
        if (opts->csize == 0) opts->csize = 0;
        if (opts->max_its==0) opts->max_its = opts->nvecs;
    } else {
        if (strcmp(task_str,"classify")==0) {
            opts->task=1;
            svm_options_read_model_file(opts);
        } else {
            if (strcmp(task_str,"translate")==0) {
                opts->task=2;
                svm_options_translate(opts);
            } else {
                fprintf(stderr,"unknown task %s in svm_options\n",task_str);
            }
        }
    }
#ifdef _OPENMP
    if (opts->nths) {
#pragma omp parallel 
        {
            omp_set_num_threads(opts->nths);
        }
    }else{
#pragma omp parallel  
        {
            if (omp_get_thread_num()==0) opts->nths = omp_get_num_threads();
        }    
    }
#endif

    fprintf(stderr,"svm options are:\n");
    svm_options_write(opts,stderr);
    program_options_free(popts);
    return opts;
}

void svm_options_free(svm_options_t* opts) {
    free(opts->vecs);
    free(opts->y);
    free(opts);
    opts = 0x0;
}

void svm_options_write(svm_options_t* opts,FILE *fp)
{
    size_t vecs_size = opts->nvecs;
    size_t kmat_size = opts->csize;
    size_t nfeat_ = opts->nfeat;
    double vsize,ksize;
    fprintf(fp,"svm options\n");
    fprintf(fp,"data file = %s\n",opts->data);
    fprintf(fp,"model file= %s\n",opts->model);
    fprintf(fp,"output    = %s\n",opts->out);
    fprintf(fp,"# nthreads  = %d\n",opts->nths);
    fprintf(fp,"task        = %d\n",opts->task);
    if (opts->task==0) {
        vsize = (sizeof(double)*vecs_size * nfeat_) / 1048576.0;
        ksize = (sizeof(double)*vecs_size * kmat_size)/1048576.0;
        fprintf(fp,"# vectors  = %d \n",opts->nvecs);
        fprintf(fp,"# features = %d \n",opts->nfeat);
        fprintf(fp,"scale kernel= %d\n",opts->scale_kernel);
        fprintf(fp,"kernel type = %d\n",opts->ktype);
        fprintf(fp,"kernel pow  = %d\n",opts->kpow);
        fprintf(fp,"kernel c1   = %le\n",opts->kc1);
        fprintf(fp,"kernel c2   = %le\n",opts->kc2);
        fprintf(fp,"cache size  = %lu\n",opts->csize);
        fprintf(fp,"eps         = %le\n",opts->eps);
        fprintf(fp,"cost        = %le\n",opts->cost);
        fprintf(fp,"maxits      = %d\n",opts->max_its);
        fprintf(stderr,"vecs size   = %le MB\n",vsize);
        fprintf(stderr,"kmat size   = %le MB\n",ksize);
    } 
    return;
}

