#include "svm_options.h"


void svm_options_translate(svm_options_t* options) {
    char *pstr = strstr(options->data,".tdo");

    fprintf(stderr,"translating\n");
    fprintf(stderr,"data file is %s\n",options->data);
    fprintf(stderr,"out file is %s\n",options->out);
    svm_options_read_data_file(options);
    if (strstr(options->data,".tdo")!=0x0) {
        /* data file is a tdo file, write out libsvm file */
        write_libsvm_file(options->out,options->nvecs,options->nfeat,options->y,options->vecs);
    } else {
        write_tdo_file(options->out,options->nvecs,options->nfeat,options->y,options->vecs);
    }
}

int find_substring(const char *str,const char *sub) {
    char * s = strstr(str,sub);
    if (s) return (s-str);
    return -1;
}

void svm_options_read_data_file(svm_options_t* opts) {
    int nt;
    int nf;
    int i;
    char * px= strstr(opts->data,".tdo");

    fprintf(stderr,"reading data file %s\n",opts->data);
    if (px) {
        read_tdo_file(opts->data,&(opts->nvecs),&(opts->nfeat),&(opts->y),&(opts->vecs));
    } else {
        read_libsvm_file(opts->data,&(opts->nvecs),&(opts->nfeat),&(opts->y),&(opts->vecs));
    }
    fprintf(stderr,"read data file %s\n",opts->data);
    nt = 0;
    nf = 0;
    for (int i=0; i<opts->nvecs ; ++i) {
        if (opts->y[i]>0.0) ++nt;
        else ++nf;
    }
    fprintf(stderr,"# true = %ld #false = %ld\n",nt,nf);
}

svm_options_t * svm_options_init(int argc,char **argv) {
    int task;
    char *end;
    char *task_str;
    svm_options_t * opts;
    program_options_t * popts = program_options_init();
    program_options_insert_options(popts,"svm.options");
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
        if (opts->max_its==0) opts->max_its = opts->nvecs;
    } else {
        if (strcmp(task_str,"classify")==0) {
            opts->task=1;
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
    if (opts->nths>0) {
        omp_set_num_threads(opts->nths);
    } else {
        omp_set_num_threads(1);
        opts->nths = 1;
    }
    {
        int nt,id,nx[0];
        #pragma omp parallel private(id,nt)
        {
            id = omp_get_thread_num();
            if ( id == 0) nx[0] = omp_get_num_threads();
        }
        fprintf(stderr,"there are %d threads\n",nx[0]);
    }
#else
    fprintf(stderr,"no openmp present\n");
    opts->nths = 1;
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

void svm_options_write(svm_options_t* opts,FILE *fp) {
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

