#include "putils_c.h"

#ifdef __cplusplus
extern "C" {
#endif

void default_error_function() { exit(EXIT_FAILURE);}

static void (*err_fun)() = &default_error_function;

void putils_set_error_function( void (*new_error_function)())
{
    err_fun = new_error_function;
}

void putils_error_function()
{
    (*err_fun)();
}

void putils_system_warning(const char *what)
{
    fprintf(stderr,"%s\n",what);
}

void putils_system_error(const char *what)
{
    fprintf(stderr,"fatal error : %s %s\n",what,strerror(errno));
    putils_error_function();
}

void putils_OpenFileError(const char *name, const char *mode)
{
    fprintf(stderr,"could not open the file %s in mode %s\n",name,mode);
    putils_error_function();
}

void * Malloc(size_t n)
{
    void * ptr = malloc(n);
    if (!ptr) {
        fprintf(stderr,"Malloc cannot allocate %lu bytes\n",n);
        putils_error_function();
    }
    return ptr;
}

void * Calloc(size_t n)
{
    void * ptr = malloc(n);
    if (!ptr) {
        fprintf(stderr,"Calloc cannot allocate %lu bytes\n",n);
        putils_error_function();
    }
    memset(ptr,0x0,n);
    return ptr;
}

void * AlignedAlloc(size_t al_size,size_t n)
{
    void * ptr;
    int e = posix_memalign(ptr,al_size,n);
    if (e == -1) {
        fprintf(stderr,"Aligned alloc failed for %lu byes - alignment = %lu\n",
            n,al_size);
        fprintf(stderr,"%s\n",strerror(errno));
        putils_error_function();
    }
    return ptr;
}

void * Resize( void * ptr, size_t old_size, size_t new_size)
{
    size_t copy_size;
    void * tmp;
    if ( old_size == 0) {
        if ( ptr ) {
            free(ptr);
        }
        ptr = Malloc(new_size);
        return ptr;
    }
    if ( new_size == 0) {
        if ( ptr ) free(ptr);
        return 0x0;
    }
    if ( old_size == new_size ) return ptr;
    copy_size = (new_size > old_size) ? old_size:new_size;
    tmp = Malloc(new_size);
    if ( ptr ) {
        memcpy(tmp,ptr,copy_size);    
        free(ptr);
    }
    ptr = tmp;
    return ptr;
}

char * Strdup(const char *str) {
    char * new_str = strdup(str);
    if ( !new_str ) putils_system_error("strdup failed");
    return new_str; 
}

FILE * Fopen(const char *file_name, const char *mode)
{
    FILE * fptr = fopen(file_name,mode);
    if (!fptr) {
        fprintf(stderr,"Fopen failed for %s in mode %s\n",file_name,mode);
        putils_error_function();
    }
    return fptr;
}

FILE * Fdopen(int file_desc, const char * mode)
{
    FILE * fptr = fdopen(file_desc,mode);
    if (!fptr) {
        fprintf(stderr,"Fdopen failed for %d in mode %s\n",file_desc,mode);
        putils_error_function();
    }
    return fptr;
}

FILE * Fmemopen(char * buffer,size_t buffer_size,const char *mode)
{
    FILE * fptr = fmemopen(buffer,buffer_size,mode);
    if (!fptr) putils_system_error("fmemopen failed");
    return fptr;
}

FILE * Popen( const char *exec_cmd, const char *mode)
{
    FILE * fptr = popen(exec_cmd,mode);
    if (!fptr) {
        putils_system_error("popen failed");
    }
    return fptr;
}

void Pclose(FILE *fp)
{
    int e = pclose(fp);
    if (e) putils_system_error("pclose failed");
}

void Fwrite(const void * data, size_t obj_size, size_t cnt, FILE *fp)
{
    ssize_t nwr = fwrite(data,obj_size,cnt,fp);
    if ( nwr != cnt) {
        putils_system_error("fwrite failed");
    }
}

void Fread(void * data, size_t obj_size, size_t cnt, FILE *fp)
{
    ssize_t nrd = fwrite(data,obj_size,cnt,fp);
    if ( nrd != cnt) {
        if ( feof(fp) ) {
            fprintf(stderr,"fread reached end of file\n");
        }
        putils_system_error("fread failed");
    }
}

ssize_t Getline(char ** data, size_t *data_size, FILE *fp)
{
    ssize_t nrd = getline(data,data_size,fp);
    if (nrd == -1) {
        if ( feof(fp) ) {
            fprintf(stderr,"getline reached end of file\n");
        }
        putils_system_error("getline failed");        
    }
    return nrd;
}

ssize_t GetDelim( char ** data, size_t * data_size, char delim, FILE *fp)
{
    ssize_t nrd = getdelim(data,data_size,delim,fp);
    if (nrd == -1) {
        if ( feof(fp) ) {
            fprintf(stderr,"getdelim reached end of file\n");
        }
        putils_system_error("getdelim failed");        
    }
    return nrd;
}

int Open( const char * filename, int mode)
{
    int fdesc = open(filename,mode);
    if (fdesc == -1) putils_system_error("open failed");
    return fdesc;
}

void Read( int fdesc, void *data, size_t n)
{
    ssize_t nrd;
    char * ptr = data;
    size_t rem = n;
    while (1) {
        nrd = read(fdesc,(void*)ptr,rem);
        if ( nrd == rem ) return;
        if ( nrd == -1) {
            putils_system_error("read failed");
        } 
        ptr += nrd;
        rem -= nrd;
    }
}

void Write( int fdesc, const void *data,size_t n)
{
    ssize_t nwr;
    const char * ptr = data;
    size_t rem = n;
    while (1) {
        nwr = write(fdesc,(const void*)ptr,rem);
        if ( nwr == rem ) return;
        if ( nwr == -1) {
            putils_system_error("read failed");
        } 
        ptr += nwr;
        rem -= nwr;
    }
}

void BlockingRead( int fdesc, void *data, size_t n)
{
    int e;
    size_t nrd;
    char * ptr = data;
    size_t rem = n;
    struct timespec twait;
    struct timespec tdurr;
    double dwait = 1./CLOCKS_PER_SEC;
    double dmax = 2.;
    double acc = 0.;
    twait.tv_sec = (time_t)floor(dwait);
    twait.tv_nsec = (long)(1.e9*(dwait-twait.tv_sec));
    while (1) {
        nrd = read(fdesc,(void*)ptr,rem);
        if ( nrd == rem ) return;
        if ( nrd == -1) {
            putils_system_error("read failed");
        } 
        ptr += nrd;
        rem -= nrd;
        e = nanosleep(&twait,&tdurr);
        if (e) {
            putils_system_warning("nanosleep interrupted by signal");
        }    
        acc += dwait;
        if ( acc > dmax) {
            putils_system_error("BlockingRead : nanosleep timed out");
        }
    }
}

pid_t Fork()
{
    pid_t p = fork();
    if ( p == -1) putils_system_error("fork failed");
    return p;
}

void Pipe( int *desc)
{
    int e = pipe(desc);
    if (e) putils_system_error("pipe failed"); 
}

void putils_adjust_file_descriptors( pid_t p, int * fds1, int * fds2, int * pfds)
{
        if ( p ) {
            // parent
            close(fds1[0]);
            close(fds2[1]);
            pfds[0] = fds2[0];
            pfds[1] = fds1[1];
        }else{
            close(fds1[1]);
            close(fds2[0]);
            pfds[0] = fds1[0];
            pfds[1] = fds2[1];
        }
}

void putils_adjust_file_descriptors_one_way( pid_t p, int * fds, int * pfds, int parent_writes)
{
    if (parent_writes) {
        if ( p ) {
            // parent
            close(fds[0]);
            pfds[0] = -1;
            pfds[1] = fds[1];
        }else{
            close(fds[1]);
            pfds[0] = fds[0];
            pfds[1] = -1;;
        }
    }else{
        if ( p ) {
            // parent
            close(fds[1]);
            pfds[0] = fds[0];
            pfds[1] = -1;
        }else{
            close(fds[0]);
            pfds[0] = -1;
            pfds[1] = fds[1];
        }    
    }
}

size_t tokenize_c_string( char * str, char *delims, char ** tokens, size_t max_ntokens)
{
    char * last;
    size_t n = 0;

    char * tok = strtok_r(str,delims,&last);
    while ( tok ) {
        tokens[n] = strdup(tok);
        ++n;
        if ( n == max_ntokens ) {
            putils_system_error("tokenize_c_string max tokens exceeded");
        }
        strtok_r(0x0,delims,&last);
    }    
    return n;
}

char * c_string_to_upper(const char *str)
{
    char * new_str;
    size_t sz = strlen(str);
    size_t i;
    new_str = Strdup(str);
    for (i=0;i<sz;++i) new_str[i] = toupper(new_str[i]);
    return new_str;    
}

char * c_string_to_lower( const char *str)
{
    char * new_str;
    size_t sz = strlen(str);
    size_t i;
    new_str = Strdup(str);
    for (i=0;i<sz;++i) new_str[i] = tolower(new_str[i]);
    return new_str;    
}

char * join_c_strings(int argc,char **argv)
{
    char * jstr;
    size_t sz;
    size_t i;
    
    if (argc == 0) return 0x0;
    if ( argv == 0x0) {
        putils_system_error("join c strings failed null string");
    }
    for (i=0;i<argc;++i) {
        if ( argv[i] == 0x0) {
            putils_system_error("join c strings failed null string");
        }
        sz += strlen(argv[i]); 
    }
    sz += argc;
    jstr = (char *)Malloc(sz);
    sz = strlen(argv[0]);
    strncpy(jstr,argv[0],sz);
    for (i=1;i<argc;++i) {
        strncat(jstr," ",1);
        sz = strlen(argv[i]);
        strncat(jstr,argv[i],sz);
    }
    return jstr;
}


