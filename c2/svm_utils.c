#include "svm_utils.h"

void * Malloc(size_t n) {
#ifdef HAVE_POSIX_MEMALIGN
    size_t aln = 64; /* 256 byte width alignment */
    void *ptr;
    posix_memalign(&ptr,aln,n);
#else
    void * ptr = malloc(n);
#endif
    if (ptr) return ptr;
    fprintf(stderr,"could not allocate %zu bytes\n",n);
    exit(EXIT_FAILURE);
}

void * Calloc(size_t n) {
    void * ptr = Malloc(n);
    memset(ptr,0x0,n);
    return ptr;
}

void * Realloc(void *ptr,size_t old_size,size_t new_size) {
    void * tmp = Malloc(new_size);
    memcpy(tmp,ptr,old_size);
    free(ptr);
    ptr = tmp;
    return ptr;
}

char * Strdup(const char *str)
{
    char * r = strdup(str);
    if ( r == 0x0) {
        if ( str != 0x0) {
            fprintf(stderr,"strdup failed for string \"%s\"\n",str);
        }else{
            fprintf(stderr,"strdup failed empty string\n");
        }
        exit(EXIT_FAILURE);
    }
    return r;
}

FILE *Fopen(const char *name,const char *mode) {
    FILE *fp = fopen(name,mode);
    if (fp) return fp;
    fprintf(stderr,"could not open the file %s in mode %s\n",name,mode);
    exit(EXIT_FAILURE);
}

void Fread(void *ptr,size_t osize,size_t cnt,FILE *fp) {
    ssize_t rd = fread(ptr,osize,cnt,fp);
    if (rd!=cnt) {
        if (feof(fp)) {
            fprintf(stderr,"end of file\n");
        } else {
            fprintf(stderr,"error %s\n",strerror(errno));
        }
        fprintf(stderr,"read only %zd of %zu elements\n",rd,cnt);
        parse_error("Fread did not read all elements");
    }
}

void Fwrite(void *ptr,size_t osize,size_t cnt,FILE *fp) {
    ssize_t wr = fwrite(ptr,osize,cnt,fp);
    if (wr!=cnt) {
        fprintf(stderr,"file error %d\n",ferror(fp));
        parse_error("Fwrite did not write all elements");
    }
}

int64_t Getline(char **line_ptr,size_t *line_size,FILE *fp) {
    ssize_t n = getline(line_ptr,line_size,fp);
    if ( n >= 0) return n;
    if (feof(fp)) return -1;
    fprintf(stderr,"error getline failed %s\n",strerror(errno));
    exit(EXIT_FAILURE);
}


char ** tokens_init() {
    int i;
    char ** tokens = (char**)Malloc(SVM_MAX_TOKENS*sizeof(char*));
    for (i=0; i<SVM_MAX_TOKENS; ++i)
        tokens[i] = (char*)Malloc(SVM_MAX_LINE_SIZE);
    return tokens;
}

void tokens_free(char **tokens) {
    int i;
    for (i=SVM_MAX_TOKENS; i>0;) {
        --i;
        free(tokens[i]);
    }
    free(tokens);
}

int explode_string(char *str,const char *delims,char **tokens) {
    char * last;
    char * ptr;
    int cnt = 0;

    for ( ptr = strtok_r(str,delims,&last); ptr; ptr=strtok_r(NULL,delims,&last)) {
        strcpy(tokens[cnt],ptr);
        ++cnt;
        if (cnt == SVM_MAX_TOKENS) {
            fprintf(stderr,"token is string exceed max tokens\n");
            exit(EXIT_FAILURE);
        }
    }
    return cnt;
}

int parse_bool(char *str) {
    size_t slen = strlen(str);

    if (slen) {
        if (str[0]=='F' || str[0]=='f' || str[0]=='0') return 0;
    } else {
        return -1;
    }
    return 1;
}

void parse_error(const char *mess) {
    fprintf(stderr,"parse error : %s\n",mess);
    exit(EXIT_FAILURE);
}

void quit_error(const char *mess) {
    fprintf(stderr,"%s\n",mess);
    exit(EXIT_FAILURE);
}


void write_tdo_file(const char *file_name,int nvecs,int nfeat,double *y,double *vecs) {
    FILE *fp;
    fp = Fopen(file_name,"w");
    Fwrite((void*)&nvecs,sizeof(int),1,fp);
    Fwrite((void*)&nfeat,sizeof(int),1,fp);
    Fwrite((void*)y,sizeof(double),nvecs,fp);
    Fwrite((void*)vecs,sizeof(double),nvecs*nfeat,fp);
    fclose(fp);
}


void write_libsvm_file(const char *file_name,int nvecs,int nfeat,double *y,double *vecs) {
    int i,j;
    int index;
    int iy;
    FILE *fp = Fopen(file_name,"w");
    for (i=0; i<nvecs; ++i) {
        iy = (int) y[i];
        if (iy>=0) iy=1;
        else iy=-1;
        fprintf(fp," %d\n",iy);
        for (j=0; j<nfeat; ++j) {
            index = j + 1;
            fprintf(fp," %d:%lg",index,vecs[j+i*nfeat]);
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
    /* free buffer and tokens */
    return;
}

void read_tdo_file(const char *file_name,int *nvecs,int *nfeat,double **y,double **vecs) {
    int nv,nf;
    FILE *fp = Fopen(file_name,"r");
    Fread((void*)&nv,sizeof(int),1,fp);
    Fread((void*)&nf,sizeof(int),1,fp);
    *y = (double*)Malloc(sizeof(double)*nv);
    *vecs = (double*)Malloc(sizeof(double)*nv*nf);
    Fread((void*)*y,sizeof(double),nv,fp);
    Fread((void*)*vecs,sizeof(double),nv*nf,fp);
    *nvecs = nv;
    *nfeat = nf;
    fclose(fp);
}

void read_libsvm_file(const char *file_name,int *nvecs,int *nfeat,double **y,double **vecs) {
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
    for (iv=0; iv<nv; ++iv) {
        nrd = Getline(&sline,&line_size,fp);
        if (nrd==-1) break;
        ntokens = explode_string(sline,delims,tokens);
        y_p[iv] = strtod(tokens[0],&end_ptr);
        for (itoken = 1; itoken < ntokens; itoken+=2) {
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
