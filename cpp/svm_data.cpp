#include "svm_data.hpp"
namespace svmpack {
void svm_data::write_libsvm_data_file(const std::string& new_name, const double& tau = 2.e-15) {
    std::ofstream out(new_name.c_str());
    if (!out) {
        throw putils::FileOpenError(__FUNCTION__,new_name.c_str(),"w");
    }
    for (int i=0; i<nvecs; ++i) {
        int iy = (int)y[i];
        out << " " << iy;
        for (int j=0; j<nfeat; ++j) {
            if (vecs[j+i*nfeat]>tau) {
                out << " " << (j+1) << ":" << vecs[j+i*nfeat];
            }
        }
        out << "\n";
    }
    out.close();
}

void svm_data::write_tdo_data_file(const std::string& new_name) {
    std::ofstream out(new_name.c_str());
    if (!out) {
        throw putils::FileOpenError(__FUNCTION__,new_name.c_str(),"w");
    }
    out.write((char*)&nvecs,sizeof(int));
    out.write((char*)&nfeat,sizeof(int));
    out.write((char*)y,sizeof(double)*nvecs);
    out.write((char*)vecs,sizeof(double)*nvecs*nfeat);
    out.close();
}

void svm_data::read_libsvm_data_file(const std::string& data)
{
    nvecs = 0;
    int max_nfeat = 0;
    nfeat=0;
    const std::string delims(" :\n");
    std::vector<std::string> tokens;
    std::string sline;
    /// training task
    std::ifstream in(data.c_str());
    if (!in) {
        throw putils::FileOpenError(__FUNCTION__,data.c_str(),"r");
    }
    while (in) {
        if (!getline(in,sline)) {
            if ( in.eof()) break;
            throw putils::ReadError(__FUNCTION__,data.c_str());
        }
        if (sline.size()==0) break;
        if (sline[0]=='#') continue;
        std::size_t ntokens= putils::tokenize_string(sline,delims,tokens);
        if (ntokens==0) {
            break;
        }
        if (ntokens%2) {
            if (ntokens!=1)
            {
                auto index = std::stoi(tokens[ntokens-2]);
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
    int ivec = 0;
    while (in) {
        if (!getline(in,sline)) {
            if ( in.eof()) break;
            throw putils::ReadError(__FUNCTION__,data.c_str());
        }
        if (sline.size()==0) break;
        if (sline[0]=='#') continue;
        auto ntokens=putils::tokenize_string(sline,delims,tokens);
        if (ntokens==0) {
            break;
        }
        if (ntokens%2) {
            y[ivec] = std::stod(tokens[0]);
            for (std::size_t j=1; j<ntokens; j+=2) {
                int index = std::stoi(tokens[j]);
                assert(index!=0);
                --index;
                double value = std::stod(tokens[j+1]);
                vecs[ ivec*nfeat + index] = value;
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

void svm_data::read_tdo_data_file(const std::string& data)
{
    /// training task
    std::ifstream in;
    in.open(data.c_str());
    if (!in) {
        throw putils::FileOpenError(__FUNCTION__,data.c_str(),"r");
    }
    in.read((char*)&nvecs,sizeof(int));
    in.read((char*)&nfeat,sizeof(int));
    vecs = new double[nvecs*nfeat];
    y = new double[nvecs];
    in.read((char *)y,nvecs*sizeof(double));
    in.read((char *)vecs,nvecs*nfeat*sizeof(double));
    in.close();
}

void svm_data::read_data_file(const std::string& data) {
    std::cerr << "reading " << data << "\n";
    if ( data.find(".tdo")==std::string::npos) {
        read_libsvm_data_file(data);
    } else {
        read_tdo_data_file(data);
    }
    std::cerr << "read in data\n";
    std::cerr << "# of features = " << nfeat << "\n";
    std::cerr << "# of vectorrs = " << nvecs << "\n";
}

void svm_data::write_data_file(const std::string& data) {
    std::cerr << "writing data file " << data << "\n";
    if ( data.find(".tdo")==std::string::npos) {
        write_libsvm_data_file(data);
    } else {
        write_tdo_data_file(data);
    }
    std::cerr << " done with write\n";
}

void svm_data::edit(std::size_t nfeat_new)
{
    if (nfeat_new == nfeat) return;
    double * vnew = new double[nfeat_new*nvecs];
    if (nfeat_new > nfeat) {
        for (int k=0;k<nvecs;++k) {
            double * vn = vnew + k * nfeat_new;
            const double * vo = vecs + k * nfeat;
            for (int j=0;j<nfeat;++j) vn[j] = vo[j];
            for (int j=nfeat;j<nfeat_new;++j) vn[j] = 0.0;
        }
    }
    else
    {
        for (int k=0;k<nvecs;++k) {
            double * vn = vnew + k * nfeat_new;
            const double * vo = vecs + k * nfeat;
            for (int j=0;j<nfeat_new;++j) vn[j] = vo[j];
        }    
    }
    nfeat = nfeat_new;
    delete [] vecs;
    vecs = vnew;
}

void svm_translate( const std::string& data_in)
{
    std::size_t p = data_in.find(".tdo");
    if ( p == std::string::npos ) {
        std::string data_out = data_in;
        data_out += ".tdo";
        svm_data data(data_in);
        data.write_data_file(data_out);
    }else{
        std::string data_out = data_in.substr(0,p);
        svm_data data(data_in);
        data.write_data_file(data_out);    
    }
}

} // end namespace