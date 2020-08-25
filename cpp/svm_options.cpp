
#include "svm_options.hpp"

inline void svm_options::write_libsvm_file(std::string& new_name) {
    int i,j;
    const double tau = 2.e-14;
    std::ofstream out;
    out.open(new_name.c_str());
    for (i=0; i<nvecs; ++i) {
        int iy = (int)y[i];
        out << " " << iy;
        for (j=0; j<nfeat; ++j) {
            if (vecs[j+i*nfeat]>tau) {
                out << " " << (j+1) << ":" << vecs[j+i*nfeat];
            }
        }
        out << "\n";
    }
    out.close();
}

inline void svm_options::write_tdo_file(std::string& new_name) {
    std::ofstream out;
    out.open(new_name.c_str());
    out.write((char*)&nvecs,sizeof(int));
    out.write((char*)&nfeat,sizeof(int));
    out.write((char*)y,sizeof(double)*nvecs);
    out.write((char*)vecs,sizeof(double)*nvecs*nfeat);
    out.close();
}

inline void svm_options::read_libsvm_data_file() noexcept
{
    nvecs = 0;
    int max_nfeat = 0;
    size_t ntokens=0;
    nfeat=0;
    size_t ivec;
    int j;
    int index;
    int sz;
    const std::string delims(" :\n");
    std::vector<std::string> tokens;
    std::string sline;
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
        ntokens=explode_string(sline,delims,tokens);
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
    in.seekg(0);
    nfeat = max_nfeat;
    vecs = new double[nvecs*nfeat];
    y = new double[nvecs];
    memset(vecs,0x0,nvecs*nfeat*sizeof(double));
    std::cerr << "nfeat = " << nfeat << " # vecs = " << nvecs << "\n";
    ivec = 0;
    while (in) {
        getline(in,sline);
        if (sline.size()==0) break;
        if (sline[0]=='#') continue;
        ntokens=explode_string(sline,delims,tokens);
        if (ntokens==0) {
            break;
        }
        if (ntokens%2) {
            y[ivec] = std::stod(tokens[0]);
            for (j=1; j<ntokens; j+=2) {
                index = std::stoi(tokens[j]);
                value = std::stod(tokens[j+1]);
                vecs[ ivec*nfeat + index - 1 ] = value;
            }
            ++ivec;
        } else {
            std::cerr << "read libsvm_data_file error on line " << (nvecs+1) << "\n";
            exit(EXIT_FAILURE);
        }
        if (in.eof()) break;
    }
    in.close();
    std::cerr << "nfeat = " << nfeat << " # vecs = " << nvecs << "\n";
    std::cerr << "read libsvm file " << data << "\n";
    return;
}

inline void svm_options::read_tdo_data_file() noexcept
{
    /// training task
    std::ifstream in;
    in.open(data.c_str());
    if (!in) {
        std::cerr << "could not open " << data << "\n";
    }
    in.read((char*)&nvecs,sizeof(int));
    in.read((char*)&nfeat,sizeof(int));
    vecs = new double[nvecs*nfeat];
    y = new double[nvecs];
    in.read((char *)y,nvecs*sizeof(double));
    in.read((char *)vecs,nvecs*nfeat*sizeof(double));
    in.close();
}

inline void svm_options::read_model_file() noexcept
{
    /// training task
    std::ifstream in;
    in.open(model.c_str());
    if (!in) {
        std::cerr << "could not open " << model << "\n";
    }
    in.read((char*)&nvecs,sizeof(int));
    in.read((char*)&nfeat,sizeof(int));
    in.read((char*)&ktype,sizeof(int));
    in.read((char*)&kpow,sizeof(int));
    in.read((char*)&kc1,sizeof(double));
    in.read((char*)&kc2,sizeof(double));
    in.read((char*)&bias,sizeof(double));
    int itmp;
    in.read((char*)&(itmp),sizeof(int));
    scale_kernel = (itmp>0) ? 1:0;
    vecs = new double[nvecs*nfeat];
    y = new double[nvecs];
    in.read((char*)y,nvecs*sizeof(double));
    in.read((char*)vecs,nvecs*nfeat*sizeof(double));
    in.close();
    std::cerr << "# model vectors = " << nvecs << "\n";
    std::cerr << "# model features= " << nfeat << "\n";
    std::cerr << "model kernel function = " << ktype << "\n";
    std::cerr << "      parameters " << ktype << " " << kc1 << " " << kc2 << "\n";
    std::cerr << "bias = " << bias << "\n";
    std::cerr << "scale kernel = " << scale_kernel << "\n";
}

inline void svm_options::read_data_file() noexcept {
    if ( data.find(".tdo")==std::string::npos) {
        read_libsvm_data_file();
    } else {
        read_tdo_data_file();
    }
}

inline void svm_options::translate() {
    std::string new_file;
    if ( data.find(".tdo")==std::string::npos) {
        read_libsvm_data_file();
        new_file = data + ".tdo";
        write_tdo_file(new_file);
    } else {
        read_tdo_data_file();
        new_file = data.substr(0,data.find(".tdo"));
        write_libsvm_file(new_file);
    }
}


svm_options::svm_options(int argc,char **argv)
{
    std::string task_str;
    putils::program_options popts;

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
    popts.insert("config","config file for options","");
    popts.insert("max_iterations","max # of training cycles","0");
    popts.parse_command_line(argc,argv);
    if (popts.has_value("config"))
        popts.parse_config_file(popts.get_value("config"));

    data = popts.get_value("data");
    model = popts.get_value("model");
    out = popts.get_value("out");
    task_str = popts.get_value("task");
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
    std::cerr << "task string = " << task_str << "\n";
    task = 0;
    if (task_str.compare("train")==0) {
        read_data_file();
        if (csize==-1) csize = nvecs/6;
        if (csize==0) csize = nvecs;
        if (max_its==0) max_its = nvecs;
        task=0;
    } else {
        if (task_str.compare("classify")==0) {
            task=1;
        } else {
            if (task_str.compare("translate")==0) {
                task=2;
                translate();
            } else {
                std::cerr << "unknown task : " << std::string(task_str) << " in svm_options\n";
            }
        }
    }
    std::cerr << "svm options are:\n";
    write_file(stderr);
}


void svm_options::write_file(FILE *fp) const noexcept {
    double vsize = double(sizeof(double)* nvecs * nfeat) / 1048576.0;
    double ksize = double(sizeof(double)* nvecs * csize)/1048576.0;
    fprintf(fp,"svm options\n");
    fprintf(fp,"data file = %s\n",data.c_str());
    fprintf(fp,"model file= %s\n",model.c_str());
    fprintf(fp,"output    = %s\n",out.c_str());
    fprintf(fp,"nthreads  = %d\n",nths);
    if (task!=1) {
    fprintf(fp,"# vectors  = %d \n",nvecs);
    fprintf(fp,"# features = %d \n",nfeat);
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
        fprintf(stderr,"kmat size   = %le MB\n",ksize);
        fprintf(stderr,"vecs size   = %le MB\n",vsize);
    } 
    }
    return;
}

