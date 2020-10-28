#include "utils.hpp"


FILE *Fopen(const char *name,const char *mode) {
    FILE *fp = fopen(name,mode);
    if (fp) return fp;
    fprintf(stderr,"could not open the file %s in mode %s\n",name,mode);
    exit(EXIT_FAILURE);
}

size_t Fread(void *ptr,size_t osize,size_t cnt,FILE *fp)
{
    size_t rd = fread(ptr,osize,cnt,fp);
    if (rd!=cnt) {
        parse_error("Fread did not read all elements");
    }
    return rd;
}

size_t Fwrite(const void *ptr,size_t osize,size_t cnt,FILE *fp)
{
    size_t wr = fwrite(ptr,osize,cnt,fp);
    if (wr!=cnt) {
        parse_error("Fwrite did not write all elements");
    }
    return wr;
}

template < typename Tp >
size_t WriteStream(const Tp * ptr,size_t cnt,std::ostream& os) noexcept
{
    os.write((const char*)ptr,sizeof(Tp)*cnt);
    if (os) return cnt;
    std::cerr << "write failed for " << cnt*sizeof(Tp) << " bytes\n";
    exit(EXIT_FAILURE);
}

template < typename Tp >
size_t ReadStream(Tp * ptr,size_t cnt,std::istream& is) noexcept
{
    is.read((char*)ptr,sizeof(Tp)*cnt);
    if (is) return cnt;
    std::cerr << "read failed for " << cnt*sizeof(Tp) << " bytes\n";
    exit(EXIT_FAILURE);
}

std::istream& Getline(std::istream& is,std::string& sline)
{
    if (getline(is,sline)) return is;
    if (is.eof()) {
        return is;
    }
    if (is.bad()) {
        std::cerr << "bad bit in geline\n";
    }else{
        if (is.fail()) {
            std::cerr << "fail bit in getline\n";
        }else{
            std::cerr << "unknown failure in getline\n";
        }
    }    
    exit(EXIT_FAILURE);
}

int explode_string(const std::string& sline,const std::string& delims,std::vector<std::string>& tokens)
{
    tokens.clear();
    size_t fn,st;

    st = sline.find_first_not_of(delims,0);
    fn = sline.find_first_of(delims,st+1);
    while (st!=std::string::npos) {
        tokens.push_back(sline.substr(st,fn-st));
        st = sline.find_first_not_of(delims,fn+1);
        fn = sline.find_first_of(delims,st+1);
    }
    return tokens.size();
}

double eval0(const double* v1,const double *v2,int nfeat) noexcept
{
    double s=0.0;
    for (int i=0; i<nfeat; ++i) s+=v1[i]*v2[i];
    return s;
}

double eval1(const double* v1,const double *v2,int nfeat)  noexcept
{
    double t;
    double s=0.0;
    for (int i=0; i<nfeat; ++i) {
        t = v1[i]-v2[i];
        s+=t*t;
    }
    return s;
}

double dpowi(double x,int i) noexcept
{
    double y,z;

    if (i<0) {
        x=1./x;
        i=-i;
    }
    switch (i) {
    case 0:
        return 1.;
    case 1:
        return x;
    case 2:
        return x*x;
    case 3:
        return x*x*x;
    case 4:
        /* naive requires 3 mults we use 2 */
        y = x * x;
        return y*y;
    case 5:
        /* naive = 4 mults, use 3 */
        y = x * x;
        return x*y*y;
    case 6:
        /* naive = 5 mults, use 4 */
        y = x * x * x;
        return y*y;
    case 7:
        /* naive use 6 mults, use 5 */
        y = x * x * x;
        return y*y*x;
    case 8:
        /* naive 7 mults vs 3 */
        y = x * x;
        y = y * y;
        return y*y;
    case 9:
        /* naive 8 mults vs 4 */
        y = x * x;
        y = y * y;
        return y*y*x;
    case 10:
        /* naive 9 mults vs 4 */
        y = x * x;
        z = y * y;
        return z*z*y;
    case 11:
        /* naive 10 mults vs 5 */
        y = x * x;
        z = y * y;
        return z*z*y*x;
    case 12:
        /* naive 11 mults vs 4 */
        y = x * x;
        z = y * y;
        return z*z*z;      
    default:
        // very unlikely someone would use kpow > 12 
        // i > 12
        y = x * x;
        z = y * y;
        z = z * z * z;
        // recursive call until i-n <= 12
        // so e.g. i=25 dpowi(x,12)*dpowi(x,12)*dpowi(x,1)
        //         i=25,i=13,i=1
        return z * dpowi(x,i-12);
    }
    return 1.e307;
}

int parse_bool(const std::string& str) noexcept
{
    size_t slen = str.size();
    if (slen) {
        if (str[0]=='F' || str[0]=='f' || str[0]=='0') return 0;
    } else {
        // return -1 to indicate an error
        std::cerr << "parse bool called on empty string!\n";
        exit(EXIT_FAILURE);
    }
    return 1;
}

void analyze(int ntp,int ntn,int nfp,int nfn) noexcept
{
    const size_t nerr = nfp + nfn;
    const size_t ncorr = ntn + ntp;
    const size_t n = nerr + ncorr;
    size_t npos = ntp + nfp;
    size_t nneg = ntn + nfn;
    size_t ntrue = ntp + nfn;
    size_t nfalse = ntn + nfp;
    double dpos = npos;
    double dneg = nneg;
    double dtrue = ntrue;
    double dfalse = nfalse;
    double dn = n;
    double dc = ntp + ntn;
    double de = nfp + nfn;   
    double acc = dc/dn;
    double sens = ntp/dtrue;
    double spec = ntn/dfalse;
    double prec = ntp/dpos;
    double threat_score = double(ntp)/double(dpos + nfn);
    double bal_acc = (spec + sens ) * 0.5;  
    double f1 = (2.*ntp)/double(2.*ntp + nfp + nfn);
    double mcc = double(ntp)*double(ntn)-double(nfn)*double(nfp);
    mcc /= sqrt(dpos*dneg*dtrue*dfalse);
    std::cerr << "# predictions          = " << n << "\n";
    std::cerr << "accuracy               = " << acc << "\n";
    std::cerr << "balanced accuracy      = " << bal_acc << "\n";
    std::cerr << "# errors               = " << (nfp+nfn) << "\n";
    std::cerr << "# true                 = " << (ntp+nfn) << "\n";
    std::cerr << "# false                = " << (ntn+nfp) << "\n";
    std::cerr << "# positive             = " << (ntp+nfp) << "\n";
    std::cerr << "# negative             = " << (ntn+nfn) << "\n";
    std::cerr << "# true-positive        = " << ntp << "\n";
    std::cerr << "# true-negative        = " << ntn << "\n";
    std::cerr << "# false-positive       = " << nfp << "\n";
    std::cerr << "# false-negative       = " << nfn << "\n";
    std::cerr << "precision(ppv)         = " << prec << "\n";
    std::cerr << "sensitivity(tpr)       = " << sens  << "\n";
    std::cerr << "specificity(tnr)       = " << spec  << "\n";
    std::cerr << "Likelihood ratio pos   = " << (sens / ( 1. - spec ))  << "\n";
    std::cerr << "Likelihood ratio neg   = " << (( 1. - sens ) / spec)  << "\n";
    std::cerr << "F measure              = " << f1  << "\n";
    std::cerr << "threat score           = " << threat_score << "\n";
    std::cerr << "Matthews's Correlation = " << mcc  << "\n";
}


void parse_error(const char *mess)
{
    fprintf(stderr,"parse error : %s\n",mess);
    exit(EXIT_FAILURE);
}

void quit_error(const char *mess) {
    fprintf(stderr,"%s\n",mess);
    exit(EXIT_FAILURE);
}

