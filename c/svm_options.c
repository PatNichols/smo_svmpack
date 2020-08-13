
#include "svm_options.h"

inline void read_libsvm_data_file(svm_options_t* opts)
{
    size_t nvecs = 0;
    int max_nfeat = 0;
    size_t ntokens=0;
    size_t nfeat=0;
    size_t ivec;
    int j;
    int index;
    int sz;
    const char *delims = " :\n";
    char *end;
    size_t buffer_size = MAX_LINE_SIZE;
    char *buffer = (char*)Malloc(MAX_LINE_SIZE);
    char **tokens = tokens_init();
    double *vecs;
    double *y;
    double value;
    int i;
    double tau=1.e-14;
    FILE *fp_out;
    /// training task
    FILE *fp = Fopen(opts->data,"r");

    while (!feof(fp)) {
        sz = getline(&buffer,&buffer_size,fp);
        if (sz==0) break;
        if (buffer[0]=='#') continue;
        ntokens = explode_string(buffer,delims,tokens);
        if (ntokens==0) {
            break;
        }
        if (ntokens%2) {
            if (ntokens!=1)
            {
                index = atoi(tokens[ntokens-2]);
                if (index > max_nfeat) max_nfeat = index;
            }
            ++nvecs;
        } else {
            parse_error("read_libsvm_data_file format error");
        }
    }
    clearerr(fp);
    rewind(fp);
    nfeat = max_nfeat;
    vecs = (double*)Calloc(sizeof(double)*nvecs*nfeat);
    y = (double*)Malloc(sizeof(double)*nvecs);
    ivec = 0;
    while (!feof(fp)) {
        sz = getline(&buffer,&buffer_size,fp);
        if (sz==0) break;
        if (buffer[0]=='#') continue;
        ntokens = explode_string(buffer,delims,tokens);
        if (ntokens==0) {
            break;
        }
        y[ivec] = (double)atoi(tokens[0]);
        for (j=1; j<ntokens; ++j) {
            index = atoi(tokens[j]);
            ++j;
            value = strtod(tokens[j],&end);
            vecs[ ivec*nfeat + index - 1 ] = value;
        }
        ++ivec;
    }
    fclose(fp);
    /* free buffer and tokens */
    tokens_free(tokens);
    free(buffer);
    /* assign variables to opts */
    opts->nfeat = nfeat;
    opts->nvecs = nvecs;
    opts->y = y;
    opts->vecs = vecs;
    fprintf(stderr,"read libsvm file # vecs = %d # of featurs = %d\n",opts->nvecs,opts->nfeat);
    return;
}

inline void read_tdo_data_file(svm_options_t* opts)
{
    double *vecs;
    double *y;
    int nvecs;
    int nfeat;
    int nrd;
    /// training task
    FILE *fp = Fopen(opts->data,"r");
    Fread(&nvecs,1,sizeof(int),fp);
    Fread(&nfeat,1,sizeof(int),fp);
    vecs = (double*)Calloc(sizeof(double)*nvecs*nfeat);
    y = (double*)Malloc(sizeof(double)*nvecs);
    Fread(y,nvecs,sizeof(double),fp);
    Fread(vecs,nvecs*nfeat,sizeof(double),fp);
    fclose(fp);
    opts->nfeat = nfeat;
    opts->nvecs = nvecs;
    opts->vecs = vecs;
    opts->y = y;
    return;
}

inline void read_model_file(svm_options_t* opts)
{
    double *y;
    double *vecs;
    int nvecs;
    int nfeat;
    size_t nrd;
    /// training task
    FILE *fp = Fopen(opts->data,"r");
    Fread(&nvecs,1,sizeof(int),fp);
    Fread(&nfeat,1,sizeof(int),fp);
    Fread(&(opts->ktype),1,sizeof(int),fp);
    Fread(&(opts->kpow),1,sizeof(int),fp);
    Fread(&(opts->kc1),1,sizeof(double),fp);
    Fread(&(opts->kc2),1,sizeof(double),fp);
    Fread(&(opts->bias),1,sizeof(double),fp);
    vecs = (double*)Calloc(sizeof(double)*nvecs*nfeat);
    y = (double*)Malloc(sizeof(double)*nvecs);
    Fread(y,nvecs,sizeof(double),fp);
    Fread(vecs,nvecs*nfeat,sizeof(double),fp);
    fclose(fp);
    opts->nfeat = nfeat;
    opts->nvecs = nvecs;
    opts->vecs = vecs;
    opts->y = y;
    return;
}

int find_substring(const char *str,const char *sub)
{
    char * s = strstr(str,sub);
    if (s) return (s-str);
    return -1;
}

inline void read_data_file(svm_options_t* opts) {
    char * name = opts->data;
    char * px= strstr(name,".tdo");
    if (px==0x0)  read_libsvm_data_file(opts);
    else read_tdo_data_file(opts);
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
        read_data_file(opts);
        if (opts->csize == -1) opts->csize = opts->nvecs/6;
        if (opts->csize == 0) opts->csize = opts->nvecs;
        if (opts->max_its==0) opts->max_its = opts->nvecs;
    } else {
        if (strcmp(task_str,"classify")==0) {
            opts->task=1;
            read_model_file(opts);
        } else {
            if (strcmp(task_str,"translate")==0) {
                opts->task=2;
                translate(opts->data);
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
    fprintf(fp,"scale kernel= %d\n",opts->scale_kernel);
    fprintf(fp,"kernel type = %d\n",opts->ktype);
    fprintf(fp,"kernel pow  = %d\n",opts->kpow);
    fprintf(fp,"kernel c1   = %le\n",opts->kc1);
    fprintf(fp,"kernel c2   = %le\n",opts->kc2);
    fprintf(fp,"task        = %d\n",opts->task);
    if (!(opts->task)) {
        fprintf(fp,"cache size  = %lu\n",opts->csize);
        fprintf(fp,"eps         = %le\n",opts->eps);
        fprintf(fp,"cost        = %le\n",opts->cost);
        fprintf(fp,"maxits      = %d\n",opts->max_its);
    } else {
        fprintf(fp,"bias        = %lf \n",opts->bias);
    }
    fprintf(stderr,"vecs size   = %le MB\n",vsize);
    fprintf(stderr,"kmat size   = %le MB\n",ksize);
    return;
}

