#include "svm_options.h"


inline void svm_options_read_tdo_file(svm_options_t* opts) {
    int nvecs,nfeat;
    FILE *fp;
    
    fprintf(stderr,"Reading Data file %s\n",opts->data);
    fp =Fopen(opts->data,"r");
    fread((void*)&nvecs,1,sizeof(int),fp);
    fread((void*)&nfeat,1,sizeof(int),fp);
    fprintf(stderr,"# vecs = %ld # feat = %ld \n",nvecs,nfeat);
    opts->nvecs =nvecs;
    opts->nfeat =nfeat;
    opts->y = (double*)Malloc(sizeof(double)*nvecs);
    opts->vecs = (double*)Malloc(sizeof(double)*nvecs*nfeat);
    fread((void*)(opts->y),nvecs,sizeof(double),fp);
    fread((void*)(opts->vecs),nvecs*nfeat,sizeof(double),fp);
    fclose(fp);
    fprintf(stderr,"# vecs = %ld # feat = %ld \n",nvecs,nfeat);
}

inline void svm_options_read_libsvm_file(svm_options_t* opts)
{
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
    fprintf(stderr,"Reading Data file %s\n",opts->data);
    fp = Fopen(opts->data,"r");
    while (!feof(fp)) {
        if (!getline(&buffer,&buffer_size,fp)) break;
        if (buffer[0]=='#') continue;
        ntokens = explode_string(buffer,delims,tokens);
        if (ntokens==0) break;
        if (ntokens%2) {
            index = ntokens-2;
            sz = atoi(tokens[index]);
            if (sz > nfeat) nfeat=(int)sz; 
        }else{
            fprintf(stderr,"Format Error reading data file %s\n",opts->data);
            exit(EXIT_FAILURE);
        }
        ++nvecs;
    }
    clearerr(fp);
    rewind(fp);
    fprintf(stderr,"# vecs = %ld # feat = %ld \n",nvecs,nfeat);
    opts->nvecs =nvecs;
    opts->nfeat =nfeat;
    opts->y = (double*)Malloc(sizeof(double)*nvecs);
    opts->vecs = (double*)Malloc(sizeof(double)*nvecs*nfeat);
    memset(opts->vecs,0x0,sizeof(double)*nvecs*nfeat);
    ivec = 0;
    while (!feof(fp)) {
        if (!getline(&buffer,&buffer_size,fp)) break;
        if (buffer[0]=='#') continue;
        ntokens = explode_string(buffer,delims,tokens);
        if (ntokens==0) break;
        (opts->y)[ivec] = strtod(tokens[0],&end);
        for (j=1;j<ntokens;j+=2) {
            index = atoi(tokens[j]);
            value = strtod(tokens[j+1],&end);
            (opts->vecs)[ivec*nfeat + index-1] = value;        
        }
        ++ivec;
    }
    fclose(fp);    
    tokens_free(tokens);
    free(buffer);
}    

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
    fread((void*)&nvecs,1,sizeof(int),fp);
    fread((void*)&nfeat,1,sizeof(int),fp);
    fread((void*)&(opts->ktype),1,sizeof(int),fp);
    fread((void*)&(opts->kpow),1,sizeof(int),fp);
    fread((void*)&(opts->kc1),1,sizeof(double),fp);
    fread((void*)&(opts->kc2),1,sizeof(double),fp);
    fread((void*)&(opts->bias),1,sizeof(double),fp);
    fread((void*)&(opts->scale_kernel),1,sizeof(int),fp);
    vecs = (double*)Malloc(sizeof(double)*nvecs*nfeat);
    y = (double*)Malloc(sizeof(double)*nvecs);
    fread((void*)y,nvecs,sizeof(double),fp);
    fread((void*)vecs,nvecs*nfeat,sizeof(double),fp);
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
    int i,j;
    char * data = options->data;
    char * new_name;
    int nvecs;
    int nfeat;
    char *pstr = strstr(options->data,".tdo");

    svm_options_read_data_file(options);
    nfeat = options->nfeat;
    nvecs = options->nvecs;
    fprintf(stderr,"translating\n");
    fprintf(stderr,"# vecs = %ld\n",options->nvecs);
    fprintf(stderr,"# feat = %ld\n",options->nfeat);
    if (!pstr) {        
        /* data is not a tdo file */
        fprintf(stderr," %s is not a tdo file\n",options->data);
        new_name = strdup(options->data);
        new_name = strcat(new_name,".tdo");
        fprintf(stderr,"new name is %s\n",new_name);
        write_tdo_file(new_name,nvecs,nfeat,options->y,options->vecs);
    }else{
        /* data is a tdo file */
        fprintf(stderr," %s is a tdo file\n",options->data);
        ptrdiff_t p = pstr-data;
        new_name = strndup(options->data,p);
        fprintf(stderr,"new name is %s\n",new_name);
        write_libsvm_file(new_name,nvecs,nfeat,options->y,options->vecs);
    }
    {
        size_t sz = sizeof(double)*(nvecs+nvecs*nfeat)+sizeof(int)*2;
        fprintf(stderr,"tdo file size is %ld\n",sz);
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
     
    if (px==0x0)  {
        fprintf(stderr,"libsvm data file\n");
        svm_options_read_libsvm_file(opts);
    }
    else {
        fprintf(stderr,"tdo data file\n");
        svm_options_read_tdo_file(opts);
    }
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
    program_options_insert(popts,"nthreads","number of threads","1");
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
        if (opts->csize == 0) opts->csize = opts->nvecs;
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

    double vsize = (sizeof(double)*vecs_size * nfeat_) / 1048576.0;
    double ksize = (sizeof(double)*vecs_size * kmat_size)/1048576.0;
    fprintf(fp,"svm options\n");
    fprintf(fp,"data file = %s\n",opts->data);
    fprintf(fp,"model file= %s\n",opts->model);
    fprintf(fp,"output    = %s\n",opts->out);
    fprintf(fp,"nthreads  = %d\n",opts->nths);
    fprintf(fp,"# vectors  = %d \n",opts->nvecs);
    fprintf(fp,"# features = %d \n",opts->nfeat);
    fprintf(fp,"# threads  = %d \n",opts->nths);
    if (opts->task!=2) {
        fprintf(fp,"scale kernel= %d\n",opts->scale_kernel);
        fprintf(fp,"kernel type = %d\n",opts->ktype);
        fprintf(fp,"kernel pow  = %d\n",opts->kpow);
        fprintf(fp,"kernel c1   = %le\n",opts->kc1);
        fprintf(fp,"kernel c2   = %le\n",opts->kc2);
    }
    fprintf(fp,"task        = %d\n",opts->task);
    if (opts->task==0) {
        fprintf(fp,"cache size  = %lu\n",opts->csize);
        fprintf(fp,"eps         = %le\n",opts->eps);
        fprintf(fp,"cost        = %le\n",opts->cost);
        fprintf(fp,"maxits      = %d\n",opts->max_its);
        fprintf(stderr,"vecs size   = %le MB\n",vsize);
        fprintf(stderr,"kmat size   = %le MB\n",ksize);
    } 
    if (opts->task==1) {
        fprintf(fp,"bias        = %lf \n",opts->bias);
        fprintf(stderr,"vecs size   = %le MB\n",vsize);
    }
    return;
}

