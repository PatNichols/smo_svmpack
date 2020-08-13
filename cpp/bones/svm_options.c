
#include "svm_options.h"

inline void svm_options::read_libsvm_data_file()
{
    size_t nvecs = 0;
    int max_nfeat = 0;
    size_t ntokens=0;
    size_t nfeat=0;
    size_t ivec;
    int j;
    int index;
    int sz;
    const std::string delims(" :\n");
    std::vector<std::string> tokens;
    std::string sline;
    double *vecs;
    double *y;
    double value;
    int i;
    double tau=1.e-14;

    /// training task
    std::ifstream in;
    in.open(data.c_str());

    while (in) {
        getline(in,sline);
        if (sline.size()==0) break;
        if (sline[0]=='#') continue;
        explode_string(sline,delims,tokens);
        ntokens = tokens.size();
        if (ntokens==0) {
            break;
        }
        if (ntokens%2) {
            if (ntokens!=1)
            {
                index = std::stoi(tokens[ntokens-2]);
                if (index > max_nfeat) max_nfeat = index;
            }
            ++nvecs;
        } else {
            std::cerr << "read libsvm_data_file error on line " << (nvecs+1) << "\n"; 
            exit(EXIT_FAILURE);        
        }
        if (in.eof()) break;
    }
    in.clear();
    in.rewind();
    nfeat = max_nfeat;
    vecs = new doubloe[nvecs*nfeat];
    y = new double[nvecs];
    memset(vecs,0x0,nvecs*nfeat*sizeof(double));
    ivec = 0;
    while (in) {
        getline(in,sline);
        if (sline.size()==0) break;
        if (sline[0]=='#') continue;
        explode_string(buffer,delims,tokens);
        ntokens = token.size();
        if (ntokens==0) {
            break;
        }
        y[ivec] = std::strtod(tokens[0]);
        for (j=1; j<ntokens; ++j) {
            index = std::stoi(tokens[j]);
            ++j;
            value = std::stod(tokens[j]);
            vecs[ ivec*nfeat + index - 1 ] = value;
        }
        ++ivec;
        if (in.eof()) break;
    }
    in.close();
    /* assign variables to opts */
    nfeat = nfeat;
    nvecs = nvecs;
    fprintf(stderr,"read libsvm file # vecs = %d # of featurs = %d\n",nvecs,nfeat);
    return;
}

inline void read_tdo_data_file(svm_options* opts)
{
    /// training task
    std::ifstream in;
    in.open(data.c_str());
    if (!in) {
        std::cerr << "could not open " << data << "\n";
    }    
    in.read((void*)&nvecs,sizeof(int));
    in.read((void*)&nfeat,sizeof(int));
    vecs = new double[nvecs*nfeat];
    y = new double[nvecs];
    in.read((void *)y,nvecs*sizeof(double));
    in.read((void *)vecs,nvecs*nfeat*sizeof(double));
    in.close();
}

inline void read_model_file(svm_options* opts)
{
    /// training task
    std::ifstream in;
    in.open(data.c_str());
    if (!in) {
        std::cerr << "could not open " << data << "\n";
    }    
    in.read((void*)&nvecs,1*sizeof(int));
    in.read((void*)&nfeat,1*sizeof(int));
    in.read((void*)&(opts->ktype),1*sizeof(int));
    in.read((void*)&(opts->kpow),1*sizeof(int));
    in.read((void*)&(opts->kc1),1*sizeof(double));
    in.read((void*)&(opts->kc2),1*sizeof(double));
    in.read((void*)&(opts->bias),1*sizeof(double));
    vecs = new double[nvecs*nfeat];
    y = new double[nvecs];
    in.read((void*)y,nvecs*sizeof(double));
    in.read((void*)vecs,nvecs*nfeat*sizeof(double));
    in.close();
}

inline void svm_options::read_data_file() {
    if ( data.find(".tdo")==std::string::npos) {
        read_libsvm_data_file(opts);
    } else {
        read_tdo_data_file(opts);
    }
}

svm_options::svm_options(int argc,char **argv)
{
    int task;
    char *end;
    char *task_str;
    program_options_t popts;

    popts.insert("data","input data file","svm.in");
    popts.insert("model","model file name","svm.model");
    popts.insert("out","output file name","svm.out");
    popts.insert("task","task to perform (train,classify,translate)","train");
    popts.insert("kernel_type","kernel_function type:\n   0)dot product\n   1)polynomial dot product\
        \n   2)radial basis function\n   3)logistic function","2");
    popts.insert("kernel_power","kernel power for type 1 kernel function","2");
    popts.insert("kernel_cof1","first parameter for kernel function","-1.0");
    popts.insert("kernel_cof2","second parameter for kernel function","0.0");
    popts.insert("cost","cost parameter for soft margin training","1.0");
    popts.insert("eps","convergence parameter for training","1.e-12");
    popts.insert("nthreads","number of threads","1");
    popts.insert("cache_size","number of rows to cache","-1");
    popts.insert("scale","scale kernel so diagonal elements are 1","true");
    popts.insert("config","config file for options",0x0);
    popts.insert("max_iterations","max # of training cycles","0");
    popts.parse_command_line(argc,argv);
    if (popts.has_value("config"))
        popts.parse_config_file(popts.get_value("config"));

        
    data = (popts.get_value("data"));
    model = (popts.get_value("model"));
    out = (popts.get_value("out"));
    task_str = (popts.get_value("task"));
    ktype = std::stoi(popts.get_value("kernel_type"));
    kpow = std::stoi(popts.get_value("kernel_power"));
    kc1 = std::stod(popts.get_value("kernel_cof1"));
    kc2 = std::stod(popts.get_value("kernel_cof2"));
    cost = std::stod(popts.get_value("cost"));
    eps = std::stod(popts.get_value("eps"));
    nths = std::stoi(popts.get_value("nthreads"));
    csize = std::stoi(popts.get_value("cache_size"));
    scale_kernel = parse_bool(popts.get_value("scale"));
    max_its = std::stoi(popts.get_value("max_iterations"));

    if (task_str.compare("train")==0) {
        task=0;
        read_data_file(opts);
        if (csize == -1) csize = nvecs/6;
        if (csize == 0) csize = nvecs;
        if (max_its==0) max_its = nvecs;
    } else {
        if (task_str.compare("classify")==0) {
            task=1;
            read_model_file(opts);
        } else {
            if (task_str.compare("translate")==0) {
                task=2;
                translate(data);
            } else {
                fprintf(stderr,"unknown task %s in svm_options\n",task_str);
            }
        }
    }
    fprintf(stderr,"svm options are:\n");
    svm_options_write(opts,std::cerr);
}

svm_options::~svm_options() {
    free(vecs);
    free(y);
}

void svm_options::write(FILE *fp)
{
    size_t vecs_size = nvecs;
    size_t kmat_size = csize;

    double vsize = (sizeof(double)*vecs_size * nfeat) / 1048576.0;
    double ksize = (sizeof(double)*vecs_size * kmat_size)/1048576.0;
    fprintf(fp,"svm options\n");
    fprintf(fp,"data file = %s\n",data.c_str());
    fprintf(fp,"model file= %s\n",model.c_str());
    fprintf(fp,"output    = %s\n",out.c_str());
    fprintf(fp,"nthreads  = %d\n",nths);
    fprintf(fp,"# vectors  = %d \n",nvecs);
    fprintf(fp,"# features = %d \n",nfeat);
    fprintf(fp,"# threads  = %d \n",nths);
    fprintf(fp,"scale kernel= %d\n",scale_kernel);
    fprintf(fp,"kernel type = %d\n",ktype);
    fprintf(fp,"kernel pow  = %d\n",kpow);
    fprintf(fp,"kernel c1   = %le\n",kc1);
    fprintf(fp,"kernel c2   = %le\n",kc2);
    fprintf(fp,"task        = %d\n",task);
    if (!(task)) {
        fprintf(fp,"cache size  = %lu\n",csize);
        fprintf(fp,"eps         = %le\n",eps);
        fprintf(fp,"cost        = %le\n",cost);
        fprintf(fp,"maxits      = %d\n",max_its);
    } else {
        fprintf(fp,"bias        = %lf \n",bias);
    }
    fprintf(stderr,"vecs size   = %le MB\n",vsize);
    fprintf(stderr,"kmat size   = %le MB\n",ksize);
    return;
}

