#include "utils.h"

//#define SVM_ALIGN_ALLOC

void * Malloc(size_t n)
{
#ifdef HAVE_POSIX_MEMALIGN
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


FILE *Fopen(const char *name,const char *mode) {
    FILE *fp = fopen(name,mode);
    if (fp) return fp;
    fprintf(stderr,"could not open the file %s in mode %s\n",name,mode);
    exit(EXIT_FAILURE);
}

size_t Fread(void *ptr,size_t osize,size_t cnt,FILE *fp)
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

size_t Fwrite(void *ptr,size_t osize,size_t cnt,FILE *fp)
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


char ** tokens_init() {
    int i;
    char ** tokens = (char**)Malloc(MAX_TOKENS*sizeof(char*));
    for (i=0; i<MAX_TOKENS; ++i) tokens[i] = (char*)Malloc(MAX_LINE_SIZE);
    return tokens;
}

void tokens_free(char **tokens)
{
    int i;
    for (i=MAX_TOKENS; i>0;) {
        --i;
        if ( tokens[i] ) free(tokens[i]);
    }
    free(tokens);
}

char * string_alloc() {
    return (char*)Malloc(MAX_LINE_SIZE);
}


int explode_string(char *str,const char *delims,char **tokens)
{
    char * last;
    char * ptr;
    int cnt = 0;

    ptr = strtok_r(str,delims,&last);
    while ( ptr )
    {
       strcpy(tokens[cnt],ptr);
       ++cnt;
       ptr = strtok_r(0x0,delims,&last);  
    } 
    return cnt;
}

int parse_bool(char *str)
{
    size_t slen = strlen(str);

    if (slen) {
        if (str[0]=='F' || str[0]=='f' || str[0]=='0') return 0;
    } else {
        return 0;
    }
    return 1;
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


void write_tdo_file(const char *file_name,int nvecs,int nfeat,
    const double *y,const double *vecs)
{
    FILE *fp;    
    fp = Fopen(file_name,"w");
    Fwrite((void*)&nvecs,sizeof(int),1,fp);
    Fwrite((void*)&nfeat,sizeof(int),1,fp);
    Fwrite((void*)y,sizeof(double),nvecs,fp);
    Fwrite((void*)vecs,sizeof(double),nvecs*nfeat,fp);
    fclose(fp);
}


void write_libsvm_file(const char *file_name,int nvecs,int nfeat,
    const double *y,const double *vecs)
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
            value = vecs[j+i*nfrat];
            if (fabs(value) > 1.e-14)
            index = j + 1;
            fprintf(fp," %d:%lg",index,value);
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
    return;
}

void read_tdo_file(const char *file_name,int *nvecs,int *nfeat,double **y,double **vecs)
{
    int nv,nf;
    double *yx;
    double *vecsx;
    FILE *fp = Fopen(file_name,"r");
    Fread((void*)&nv,sizeof(int),1,fp);
    Fread((void*)&nf,sizeof(int),1,fp);
    yx = (double*)Malloc(sizeof(double)*nv);
    vecsx = (double*)Malloc(sizeof(double)*nv*nf);
    Fread((void*)yx,sizeof(double),nv,fp);
    Fread((void*)vecsx,sizeof(double),nv*nf,fp);
    *nvecs = nv;
    *nfeat = nf;
    *vecs = vecsx;
    *y = yx;
    fclose(fp);
}

void read_libsvm_file(const char *file_name,int *nvecs,int *nfeat,double **y,double **vecs)
{
    ssize_t nrd;
    size_t line_size = MAX_LINE_SIZE;
    int ntokens;
    int i,j;
    long index;
    int itoken;
    int iy;
    int iv;
    int nv = 0;
    long nf = 0;
    double *vec_p;
    double *y_p;
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
    vec_p = *vecs;
    y_p = *y;
    memset(vec_p,0x0,sizeof(double)*nv*nf);
    for (iv=0;iv<nv;++iv) {
        nrd = Getline(&sline,&line_size,fp);
        if (nrd==-1) break;
        ntokens = explode_string(sline,delims,tokens);
        y_p[iv] = strtod(tokens[0],&end_ptr); 
        for (itoken = 1;itoken < ntokens;itoken+=2) {
            index = strtoull(tokens[itoken],&end_ptr,10);
            value = strtod(tokens[itoken+1],&end_ptr);
            vec_p[iv*nf+index-1] = value; 
        }
    }
    clearerr(fp);
    fclose(fp);    
    tokens_free(tokens);
    free(sline);
}

