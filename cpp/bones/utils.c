#include "utils.h"

void translate(char * file_in)
{
    char *tdo_file;
    char *svm_file;
    if (strstr(file_in,".tdo")==0) {
        // this should be a libsvm file
        tdo_file = strcat(file_in,".tdo");
        translate_to_tdo(file_in,tdo_file);
    } else {
        svm_file = strdup(file_in);
        int i=0;
        for (;; ++i) {
            if (svm_file[i]=='.') {
                svm_file[i]=0x0;
                break;
            }
        }
        translate_to_libsvm(svm_file,file_in);
    }
}


void translate_to_tdo(const char *svm_file,const char *tdo_file)
{
    size_t buffer_size = MAX_LINE_SIZE;
    char * buffer = (char*)Malloc(MAX_LINE_SIZE);
    char ** tokens = tokens_init();
    const char *delims=" :=\n";
    FILE *fp_in;
    FILE *fp_out;
    int nfeat,max_nfeat,j,ivec,rd,ntokens;
    size_t index;
    double value;
    int nvecs = 0;
    max_nfeat=0;
    double *vecs;
    double *y;
    char *end;

    fp_in = Fopen(svm_file,"r");
    while (!feof(fp_in)) {
        rd = getline(&buffer,&buffer_size,fp_in);
        if (rd==0) continue;
        if (buffer[0]=='#') continue;
        ntokens = explode_string(buffer,delims,tokens);
        if (ntokens==0) continue;
        nfeat = atoi(tokens[ntokens-2]);
        if (nfeat > max_nfeat) max_nfeat = nfeat;
        ++nvecs;
    }
    nfeat = max_nfeat;
    clearerr(fp_in);
    rewind(fp_in);

    vecs = (double*)Calloc(sizeof(double)*nvecs*nfeat);
    y = (double*)Malloc(sizeof(double)*nvecs);

    for (ivec = 0; ivec < nvecs; ++ivec) {
        rd = getline(&buffer,&buffer_size,fp_in);
        if (rd==0) continue;
        if (buffer[0]=='#') continue;
        ntokens = explode_string(buffer,delims,tokens);
        if (ntokens==0) continue;
        y[ivec] = strtod(tokens[0],&end);
        for (j=1; j<ntokens; ++j) {
            index = atoi(tokens[j]);
            ++j;
            value = strtod(tokens[j],&end);
            vecs[index+ivec*nfeat] = value;
        }
    }
    fclose(fp_in);
    fp_out = Fopen(tdo_file,"w");
    fwrite(&nvecs,1,sizeof(int),fp_out);
    fwrite(&nfeat,1,sizeof(int),fp_out);
    fwrite(y,nvecs,sizeof(double),fp_out);
    fwrite(vecs,nvecs*nfeat,sizeof(double),fp_out);
    fclose(fp_out);
    free(y);
    free(vecs);
    tokens_free(tokens);
    free(buffer);
}

void translate_to_libsvm(const char *svm_file,const char *tdo_file)
{
    FILE *fp_in;
    FILE *fp_out;

    int nvecs,nfeat;
    double *y;
    double *vecs;
    int j,ivec;
    const double *vp;
    const double tau = 1.e-12;
    const char *delims=" :=\n";

    fp_in = Fopen(tdo_file,"r");
    fread(&nvecs,1,sizeof(int),fp_in);
    fread(&nfeat,1,sizeof(int),fp_in);
    y = (double*)Malloc(sizeof(double)*nvecs);
    vecs = (double*)Malloc(sizeof(double)*nvecs*nfeat);
    fread(y,nvecs,sizeof(double),fp_in);
    fread(vecs,nvecs*nfeat,sizeof(double),fp_in);
    fclose(fp_in);
    fp_out = Fopen(svm_file,"w");
    for (ivec=0; ivec<nvecs; ++ivec) {
        fprintf(fp_out," %2d",(int)rint(y[ivec]));
        vp = vecs + ivec*nfeat;
        for (j=0; j<nfeat; ++j) {
            if (fabs(vp[j]) > tau) {
                fprintf(fp_out," %d:%lf",(j+1),vp[j]);
            }
        }
        fprintf(fp_out,"\n");
    }
    fclose(fp_out);
    free(y);
    free(vecs);
}

#define SVM_ALIGN_ALLOC

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

inline int64_t find_str(char *str,const char *sub_str) {
    char *p  = strstr(str,sub_str);
    if (p) {
        return (p-str);
    }
    return -1;
}

inline int64_t find_first_of(char *str,const char *delims,int p)
{
    int64_t k=0;
    char ch;
    char ich;
    char *p0;
    const char *p1;
    p0 = str;
    while (str[k]!=0x0) {
        ch = str[k];
        p1 = delims;
        while (*p1!=0x0) {
            if (ch==*p1) return k;
            ++p1;
        }
        ++k;
    }
    return -1;
}

inline int64_t find_first_not_of(const char *str,const char *delims,int p)
{
    int64_t i,j;
    int is_inc;
    char ch;
    int slen = strlen(str);
    int dlen = strlen(delims);

    for (i=p; i<slen; ++i) {
        ch = str[i+p];
        is_inc = 0;
        for (j=0; j<dlen; ++j) {
            if (ch==delims[j]) {
                is_inc = 1;
                break;
            }
        }
        if (is_inc==0) return (i+p);
    }
    return -1;
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
    double mcd = sqrt ( dnp ) * dnn * dnt * dnf;
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

