#ifndef SVM_CLASSIFY_HPP
#define SVM_CLASSIFY_HPP

#include <cstdio>
#include <cstdlib>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "utils.hpp"
#include "svm_options.hpp"
#include "stopwatch.hpp"

struct svm_model {
    double *yalf;
    double *vecs;
    double c1,c2,bias;
    int nfeat,nvecs;
    int ktype,kpower,kscale,owns;

    svm_model() = default;

    svm_model(const std::string& mod)
    {
        std::ifstream in(mod.c_str());
        if (!in) {
            std::cerr << "Could not open the model file " << mod << "\n";
            exit(EXIT_FAILURE);
        }
        in.read((char*)&nvecs,sizeof(int));
        in.read((char*)&nfeat,sizeof(int));
        in.read((char*)&ktype,sizeof(int));
        in.read((char*)&kpower,sizeof(int));
        in.read((char*)&c1,sizeof(double));
        in.read((char*)&c2,sizeof(double));
        in.read((char*)&bias,sizeof(double));
        in.read((char*)&kscale,sizeof(int));
        yalf = new double[nvecs];
        vecs = new double[nvecs*nfeat];
        in.read((char*)yalf,sizeof(double)*nvecs);
        in.read((char*)vecs,sizeof(double)*nvecs*nfeat);
        in.close();
        if (ktype>3 || ktype<0) {
            std::cerr << "unknown kernel type " << ktype << "\n";
            exit(EXIT_FAILURE);
        }
        if (ktype==2) c1 = -c1;
        std::cerr << "svm model file read " << mod << "\n";
        std::cerr << "nvecs_model = " << nvecs << "\n";
        std::cerr << "nfeat_model = " << nfeat << "\n";
        int nt = 0;
        int nf = 0;
        for (int i=0; i<nvecs; ++i) {
            if (yalf[i]>=0) ++nt;
            else ++nf;
        }
        std::cerr << "# true = " << nt << " #false = " << nf << "\n";
        owns = true;
    }

    explicit svm_model(svm_model& model):
        yalf(model.yalf),vecs(model.vecs),
        c1(model.c1),c2(model.c2),bias(model.bias),
        nvecs(model.nvecs),nfeat(model.nfeat),
        ktype(model.ktype),kpower(model.kpower),kscale(model.kscale),owns(0) {}

    ~svm_model() {
        if (owns) {
            delete [] vecs;
            delete [] yalf;
        }
    }

    svm_model& operator=(svm_model& model) {
        yalf = model.yalf;
        vecs = model.vecs;
        c1=model.c1;
        c2=model.c2;
        bias = model.bias;
        nvecs = model.nvecs;
        nfeat = model.nfeat;
        ktype = model.ktype;
        kpower = model.kpower;
        kscale = model.kscale;
        owns = 0;
        return *this;
    }

    constexpr double powi(double x,int m) const noexcept
    {
        double y{1.0},z{1.0};
        if (m==0) return y;
        if (m<0) {
            m = -m;
            x = 1./x;
        }
        switch (m) {
        case 0:
            return 1.;
        case 1:
            return x;
        case 2:
            return x*x;
        case 3:
            return x*x*x;
        case 4:
            y = x * x;
            return y*y;
        case 5: // 3 mult vs 4
            y = x * x;
            return y*y*x;
        case 6: // 3 mult vs 5
            y = x * x;
            return y*y*y;
        case 7: // 4 mult vs 6
            y = x * x;
            return y*y*y*x;
        case 8: // 4 mult vs 7
            y = x * x;
            y = x * x;
            return y*y;
        case 9: // 4 mult vs 8
            y = x * x;
            y = x * x;
            return y*y*x;
        case 10:  // 4 mult vs 9
            y = x * x;
            z = y * y;
            return z*z*y;
        case 11:
            y = x * x;
            z = y * y;
            return z*z*y*x;
        case 12:
            y = x * x;
            y = x * x;
            return y*y*y;
        default:
            y = 1.0;
            for (; m>0;) {
                --m;
                y *= x;
            }
            return y;
        }
    }
    constexpr double eval_dot(const double *v1,const double *v2) const noexcept {
        double s{0.0};
        for (int i=0; i<nfeat; ++i) s+=v1[i]*v2[i];
        return s;
    };
    constexpr double eval_diff(const double *v1,const double *v2) const noexcept {
        double s{0.0};
        for (int i=0; i<nfeat; ++i) {
            double t =v1[i]-v2[i];
            s+=t*t;
        }
        return s;
    };
    double scale_factor(const double *vp) const noexcept {
        if (!kscale || ktype==2) return 1.0;
        const double tau = 1.e-15;
        double s;
        switch (ktype) {
        case 0:
            s = eval_dot(vp,vp);
            s = (s>tau) ? 1./sqrt(s):1.0;
            return s;
        case 1:
            s = powi(c1*eval_dot(vp,vp)+c2,kpower);
            s = (s>tau) ? 1./sqrt(s):1.0;
            return s;
        case 2:
            return 1.0;
        case 3:
            s = tanh(c1*eval_dot(vp,vp)+c2);
            s = (s>tau) ? 1./sqrt(s):1.0;
            return s;
        }
        return 0.0;
    }
    double gensum(const double *vp) const noexcept {
        int i;
        double scal_p = scale_factor(vp);
        double sum = 0.0;
        switch (ktype) {
        case 0:
            for (i=0; i<nvecs; ++i) {
                sum += eval_dot(vp,vecs+i*nfeat) * yalf[i];
            }
            return sum*scal_p-bias;
        case 1:
            for (i=0; i<nvecs; ++i) {
                sum += powi(c1*eval_dot(vp,vecs+i*nfeat)+c2,kpower) * yalf[i];
            }
            return sum*scal_p-bias;
        case 2:
            for (i=0; i<nvecs; ++i) {
                sum += exp( c1 * eval_diff(vp,vecs+i*nfeat) )* yalf[i];
            }
            return sum-bias;
        case 3:
            for (i=0; i<nvecs; ++i) {
                sum += tanh(c1*eval_dot(vp,vecs+i*nfeat)+c2) * yalf[i];
            }
            return sum*scal_p-bias;
        }
        return -1.e100;
    }
};


struct svm_data {
    double *y;
    double *vecs;
    int nfeat,nvecs;

    svm_data(std::string& str) {
        if (str.find(".tdo")!=std::string::npos) {
            std::ifstream in(str.c_str());
            if (!in) {
                std::cerr << "could not open the data file " << str << "\n";
            }
            in.read((char*)&nvecs,sizeof(int));
            in.read((char*)&nfeat,sizeof(int));
            y = new double[nvecs];
            vecs = new double[nvecs*nfeat];
            in.read((char*)y,sizeof(double)*nvecs);
            in.read((char*)vecs,sizeof(double)*nvecs*nfeat);
            in.close();
        } else {
            read_libsvm(str);
        }
    }

    ~svm_data() {
        delete [] vecs;
        delete [] y;
    }

    void read_libsvm(std::string& data) {
        nvecs = 0;
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
        if (!in) {
            std::cerr << "could not open the data file " << data << "\n";
        }
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
                    if (index > nfeat) nfeat = index;
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
        std::cerr << "svm data\n";
        std::cerr << "nfeat = " << nfeat << " # vecs = " << nvecs << "\n";
        std::cerr << "read libsvm file " << data << "\n";
    }
};

void svm_classify(svm_options& opts);

#endif