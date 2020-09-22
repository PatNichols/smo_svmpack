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

inline size_t Fread(void *ptr,size_t osize,size_t cnt,FILE *fp)
{
    ssize_t rd = fread(ptr,osize,cnt,fp);
    if (rd!=cnt) {
        if (feof(fp)) {
            fprintf(stderr,"end of file\n");
        }else{
            fprintf(stderr,"error %d\n",errno);
        }
        fprintf(stderr,"read only %d of %d elements\n",rd,cnt);
        parse_error("Fread did not read all elements");
    }
    return rd;
}

inline size_t Fwrite(void *ptr,size_t osize,size_t cnt,FILE *fp)
{
    ssize_t wr = fwrite(ptr,osize,cnt,fp);
    if (wr!=cnt) {
        parse_error("Fwrite did not write all elements");
    }
    return wr;
}

ssize_t Getline(char **line_ptr,size_t *line_size,FILE *fp)
{
   int p;
   ssize_t n = getline(line_ptr,line_size,fp); 
   if (n) return n;
   if (feof(fp)) return -1;
   p = errno;
   if (p==EINVAL) {
       fprintf(stderr,"Error in Getline line_ptr or line_size is null\n");
   }else{
       if (p==EOVERFLOW) {
           fprintf(stderr,"Error in Getline overflow in read\n");
       }else{
           fprintf(stderr,"Error in Getline %d\n",p);
       }
   }
   exit(EXIT_FAILURE);
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


void write_tdo_file(const char *file_name,int nvecs,int nfeat,double *y,double *vecs)
{
    FILE *fp;    
    fp = Fopen(file_name,"w");
    Fwrite((void*)&nvecs,sizeof(int),1,fp);
    Fwrite((void*)&nfeat,sizeof(int),1,fp);
    Fwrite((void*)y,sizeof(double),nvecs,fp);
    Fwrite((void*)vecs,sizeof(double),nvecs*nfeat,fp);
    fclose(fp);
}


void write_libsvm_file(const char *file_name,int nvecs,int nfeat,double *y,double *vecs)
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

void read_tdo_file(const char *file_name,int *nvecs,int *nfeat,double **y,double **vecs)
{
    int nv,nf;
    FILE *fp = Fopen(file_name,"r");
    fprintf(stderr,"opened file\n");
    Fread((void*)&nv,sizeof(int),1,fp);
    fprintf(stderr,"read nv = %d\n",nv);
    Fread((void*)&nf,sizeof(int),1,fp);
    fprintf(stderr,"read nv = %d nf = %d\n",nv,nf);
    *y = (double*)Malloc(sizeof(double)*nv);
    *vecs = (double*)Malloc(sizeof(double)*nv*nf);
    fprintf(stderr,"reading vecs and y nv = %d nf = %d\n",nv,nf);
    Fread((void*)*y,sizeof(double),nv,fp);
    fprintf(stderr,"read y\n");
    Fread((void*)*vecs,sizeof(double),nv*nf,fp);
    fprintf(stderr,"read vecs\n");
    *nvecs = nv;
    *nfeat = nf;
    fclose(fp);
}

void read_libsvm_file(const char *file_name,int *nvecs,int *nfeat,double **y,double **vecs)
{
    ssize_t nrd;
    size_t line_size = 512;
    int ntokens;
    int i,j;
    long index;
    int itoken;
    int iy;
    int iv;
    int nv = 0;
    long nf = 0;
    char *end_ptr = 0x0;
    char **tokens = tokens_init();
    char *sline = (char*)Malloc(line_size);
    const char * delims = " :\n";
    double value;
    FILE *fp = Fopen(file_name,"r");
    while (1) {
        nrd = Getline(&sline,&line_size,fp);
        if (nrd==-1) break;
        ntokens = explode_string(sline,delims,tokens);
        if (ntokens==0) break;
        ++nv;
        if ((ntokens%2)==0) {
            fprintf(stderr,"error in format for libsvm file line %d\n",nv);
            exit(EXIT_FAILURE);
        }
        if (ntokens > 2) {
            index = strtoull(tokens[ntokens-2],&end_ptr,10);        
            nf = (nf >= index) ? nf:index;
        }
    }
    clearerr(fp);
    fseek(fp,0L,SEEK_SET);
    *nvecs = nv;
    *nfeat = (int)nf;
    *y = (double*)Malloc(nv*sizeof(double));
    *vecs = (double*)Malloc(nv*nf*sizeof(double));
    memset(*vecs,0x0,sizeof(double)*nv*nf);
    for (iv=0;iv<nv;++iv) {
        nrd = Getline(&sline,&line_size,fp);
        if (nrd==-1) break;
        ntokens = explode_string(sline,delims,tokens);
        *y[iv] = strtoull(tokens[0],&end_ptr,10); 
        for (itoken = 1;itoken < ntokens;itoken+=2) {
            index = strtoul(tokens[itoken],&end_ptr,10);
            value = strtod(tokens[itoken+1],&end_ptr);
            *vecs[iv*nv+index-1] = value; 
        }
    }
    clearerr(fp);
    fclose(fp);    
    tokens_free(tokens);
    free(sline);
}
