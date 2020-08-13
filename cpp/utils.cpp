#include "utils.hpp"


FILE *Fopen(const char *name,const char *mode) {
    FILE *fp = fopen(name,mode);
    if (fp) return fp;
    fprintf(stderr,"could not open the file %s in mode %s\n",name,mode);
    exit(EXIT_FAILURE);
}

size_t Fread(void *ptr,size_t cnt,size_t osize,FILE *fp)
{
    size_t rd = fread(ptr,cnt,osize,fp);
    if (rd!=cnt) {
        parse_error("Fread did not read all elements");
    }
    return rd;
}

size_t Fwrite(const void *ptr,size_t cnt,size_t osize,FILE *fp)
{
    size_t wr = fwrite(ptr,cnt,osize,fp);
    if (wr!=cnt) {
        parse_error("Fwrite did not write all elements");
    }
    return wr;
}

int explode_string(const std::string& sline,const std::string& delims,std::vector<std::string>& tokens)
{
    tokens.clear();
    int cnt = 0;
    size_t fn,st;

    st = sline.find_first_not_of(delims,0);
    fn = sline.find_first_of(delims,st);
    while (st!=std::string::npos) {
        tokens.push_back(sline.substr(st,fn-st));
        ++cnt;
        st = sline.find_first_not_of(delims,fn+1);
        fn = sline.find_first_of(delims,st+1);
    }
    return cnt;
}

double eval0(const double* v1,const double *v2,int nfeat)
{
    double s=0.0;
    for (int i=0; i<nfeat; ++i) s+=v1[i]*v2[i];
    return s;
}

double eval1(const double* v1,const double *v2,int nfeat)
{
    double t;
    double s=0.0;
    for (int i=0; i<nfeat; ++i) {
        t = v1[i]-v2[i];
        s+=t*t;
    }
    return s;
}

double dpowi(double x,int i)
{
    double y;

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
        y = x * x;
        y = y * y;
        return y*y;
    case 9:
        y = x * x;
        y = y * y;
        return y*y*x;
    default:
        return pow(x,i);
    }
}

int parse_bool(const std::string& str)
{
    size_t slen = str.size();

    if (slen) {
        if (str[0]=='F' || str[0]=='f' || str[0]=='0') return 0;
    } else {
        return 0;
    }
    return 1;
}

void analyze(int ntp,int nfp,int ntn,int nfn)
{
    const size_t nerr = nfp + ntn;
    const size_t ncorr = nfn + ntp;
    const size_t n = nerr + ncorr;
    size_t np = nfp + ntp;
    size_t nn = nfn + ntn;
    size_t nt = ntp + ntn;
    size_t nf = nfn + nfp;
    double dntp = ntp;
    double dnfn = nfn;
    double dntn = ntn;
    double dnfp = nfp;
    double dnf = nf;
    double dnt = nt;
    double dnp = np;
    double dnn = nn;
    double dn = n;

    double sens = ( dntp ) / (dnt );
    double spec = ( dnfn ) / (dnf );
    double prec = ( dntp ) / (dnp );
    double recall = ( dntp ) / (dnt );
    double f = 2. * ( prec * recall ) / ( prec + recall );
    double mcn = (dntp * dnfn ) - (dnfp * dntn );
    double mcd = sqrt ( dnp * dnn * dnt * dnf );
    double acc = ((double)ncorr)/dn;
    std::cerr << "# predictions          = " << n << "\n";
    std::cerr << "# true                 = " << nt << "\n";
    std::cerr << "# false                = " << nf << "\n";
    std::cerr << "# positive             = " << np << "\n";
    std::cerr << "# negative             = " << nn << "\n";
    std::cerr << "# true-positive        = " << ntp << "\n";
    std::cerr << "# true-negative        = " << ntn << "\n";
    std::cerr << "# false-positive       = " << nfp << "\n";
    std::cerr << "# false-negative       = " << nfn << "\n";
    std::cerr << "accuracy               = " << acc << "\n";
    std::cerr << "# errors               = " << nerr << "\n";
    std::cerr << "precision              = " << prec << "\n";
    std::cerr << "recall                 = " << recall << "\n";
    std::cerr << "sensitivity(true acc)  = " << (dntp / dnt)  << "\n";
    std::cerr << "specificity(false acc) = " << (dnfn / dnf)  << "\n";
    std::cerr << "positive pred rate     = " << (dntp / dnp)  << "\n";
    std::cerr << "negative pred rate     = " << (dnfn / dnn) << "\n";
    std::cerr << "False Positive rate    = " << (dnfp  / dnf)  << "\n";
    std::cerr << "False Negative rate    = " << (dntn  / dnt)  << "\n";
    std::cerr << "Likelihood ratio pos   = " << (sens / ( 1. - spec ))  << "\n";
    std::cerr << "Likelihood ratio neg   = " << (( 1. - sens ) / spec)  << "\n";
    std::cerr << "F measure              = " << f  << "\n";
    if ( mcd > 1.e-14 ) {
        mcn = mcn / mcd;
    } else {
        mcn = 1;
    }
    std::cerr <<  "Matthews's Correlation = " << mcn  << "\n";
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

