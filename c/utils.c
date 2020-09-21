#include "utils.h"

//#define SVM_ALIGN_ALLOC

void * Malloc(size_t n)
{
#ifdef SVM_ALIGN_ALLOC
    size_t aln = 64; /* 256 byte width alignment */
    void *ptr;
    posix_memalign(&ptr,aln,n);
#else
    void * ptr = malloc(n);
#endif
    if (ptr) return ptr;
    fprintf(stderr,"could not allocate %uLL bytes\n",n);
    exit(EXIT_FAILURE);
}

void * Calloc(size_t n)
{
    void * ptr = Malloc(n);
    memset(ptr,0x0,n);
    return ptr;
}

void * Realloc(void *ptr,size_t old_size,size_t new_size)
{
    void * tmp = Malloc(new_size);
    memcpy(tmp,ptr,old_size);
    free(ptr);
    ptr = tmp;
    return ptr;
}


inline FILE *Fopen(const char *name,const char *mode) {
    FILE *fp = fopen(name,mode);
    if (fp) return fp;
    fprintf(stderr,"could not open the file %s in mode %s\n",name,mode);
    exit(EXIT_FAILURE);
}

inline size_t Fread(void *ptr,size_t cnt,size_t osize,FILE *fp)
{
    size_t rd = fread(ptr,cnt,osize,fp);
    if (rd!=cnt) {
        parse_error("Fread did not read all elements");
    }
    return rd;
}

inline size_t Fwrite(void *ptr,size_t cnt,size_t osize,FILE *fp)
{
    size_t wr = fwrite(ptr,cnt,osize,fp);
    if (wr!=cnt) {
        parse_error("Fwrite did not write all elements");
    }
    return wr;
}


inline char ** tokens_init() {
    int i;
    char ** tokens = (char**)Malloc(MAX_TOKENS*sizeof(char*));
    for (i=0; i<MAX_TOKENS; ++i) tokens[i] = (char*)Malloc(MAX_LINE_SIZE);
    return tokens;
}

inline void tokens_free(char **tokens)
{
    int i;
    for (i=MAX_TOKENS; i>0;) {
        --i;
        free(tokens[i]);
    }
    free(tokens);
}

inline int explode_string(char *str,const char *delims,char **tokens)
{
    char * last;
    char * ptr;
    int cnt = 0;

    for ( ptr = strtok_r(str,delims,&last); ptr; ptr=strtok_r(NULL,delims,&last)) {
        strcpy(tokens[cnt],ptr);
        ++cnt;
    }
    return cnt;
}

inline double eval0(const double* v1,const double *v2,int nfeat)
{
    double s=0.0;
    for (int i=0; i<nfeat; ++i) s+=v1[i]*v2[i];
    return s;
}

inline double eval1(const double* v1,const double *v2,int nfeat)
{
    double t;
    double s=0.0;
    for (int i=0; i<nfeat; ++i) {
        t = v1[i]-v2[i];
        s+=t*t;
    }
    return s;
}

inline double dpowi(double x,int i)
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

inline int parse_bool(char *str)
{
    size_t slen = strlen(str);

    if (slen) {
        if (str[0]=='F' || str[0]=='f' || str[0]=='0') return 0;
    } else {
        return 0;
    }
    return 1;
}

inline void analyze(int ntp,int ntn,int nfp,int nfn)
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
    fprintf(stderr,"# predictions          = %d\n", n);
    fprintf(stderr,"# true                 = %d\n", nt);
    fprintf(stderr,"# false                = %d\n", nf);
    fprintf(stderr,"# positive             = %d\n", np);
    fprintf(stderr,"# negative             = %d\n", nn);
    fprintf(stderr,"# true-positive        = %d\n", ntp);
    fprintf(stderr,"# true-negative        = %d\n", ntn);
    fprintf(stderr,"# false-positive       = %d\n", nfp);
    fprintf(stderr,"# false-negative       = %d\n", nfn);
    fprintf(stderr,"accuracy               = %le\n", acc);
    fprintf(stderr,"# errors               = %d\n", nerr);
    fprintf(stderr,"precision              = %le\n", prec);
    fprintf(stderr,"recall                 = %le\n", recall);
    fprintf(stderr,"sensitivity(true acc)  = %le\n", dntp / dnt );
    fprintf(stderr,"specificity(false acc) = %le\n", dnfn / dnf );
    fprintf(stderr,"positive pred rate     = %le\n", dntp / dnp );
    fprintf(stderr,"negative pred rate     = %le\n", dnfn / dnn );
    fprintf(stderr,"False Positive rate    = %le\n", dnfp  / dnf );
    fprintf(stderr,"False Negative rate    = %le\n", dntn  / dnt );
    fprintf(stderr,"Likelihood ratio pos   = %le\n", sens / ( 1. - spec ) );
    fprintf(stderr,"Likelihood ratio neg   = %le\n", ( 1. - sens ) / spec );
    fprintf(stderr,"F measure              = %le\n", f );
    if ( mcd > 1.e-14 ) {
        mcn = mcn / mcd;
    } else {
        mcn = 1;
    }
    fprintf(stderr, "Matthews's Correlation = %le\n", mcn );
}


inline void parse_error(const char *mess)
{
    fprintf(stderr,"parse error : %s\n",mess);
    exit(EXIT_FAILURE);
}

inline void quit_error(const char *mess) {
    fprintf(stderr,"%s\n",mess);
    exit(EXIT_FAILURE);
}


void write_tdo_file(char *file_name,
int nvecs,int nfeat,double *y,double *vecs)
{
    FILE *fp;    
    fp = Fopen(file_name,"w");
    fwrite((void*)&nvecs,1,sizeof(int),fp);
    fwrite((void*)&nfeat,1,sizeof(int),fp);
    fwrite((void*)y,nvecs,sizeof(double),fp);
    fwrite((void*)vecs,nvecs*nfeat,sizeof(double),fp);
    fclose(fp);
}

void write_libsvm_file(char *file_name,
int nvecs,int nfeat,double *y,double *vecs)
{	
    int i,j;
    int index;
    int iy;
    FILE *fp = Fopen(file_name,"w");
    for (i=0;i<nvecs;++i) {
        iy = (int) y[i];
        if (iy>=0) iy=1;
        else iy=-1;
        fprintf(fp," %d\n",iy);
        for (j=0;j<nfeat;++j) {
            index = j + 1;
            fprintf(fp," %d:%lg",index,vecs[j+i*nfeat]);
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
    /* free buffer and tokens */
    return;
}

