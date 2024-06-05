#pragma once
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <errno.h>
#include <math.h>
#include <ctype.h>
#include <fcntl.h>

#ifdef __cplusplus
extern "C" {
#endif

void putils_set_error_function( void (*new_error_function)());
void putils_error_function();
void putils_system_warning(const char *what);
void putils_system_error(const char *what);
void putils_OpenFileError(const char * filename, const char *mode);
void * Malloc(size_t n);
void * Calloc(size_t n);
void * AlignedAlloc(size_t al_size,size_t n);
void * Resize( void * ptr, size_t old_size, size_t new_size);
char * Strdup(const char *str);
FILE * Fopen(const char *file_name, const char *mode);
FILE * Fdopen(int file_desc, const char * mode);
FILE * Fmemopen(char * buffer,size_t buffer_size,const char *mode);
FILE * Popen( const char *exec_cmd, const char *mode);
void Pclose(FILE *fp);
void Fwrite(const void * data, size_t obj_size, size_t cnt, FILE *fp);
void Fread(void * data, size_t obj_size, size_t cnt, FILE *fp);
ssize_t Getline(char ** data, size_t *data_size, FILE *fp);
ssize_t GetDelim( char ** data, size_t * data_size, char delim, FILE *fp);
int Open( const char * filename, int mode);
void Read( int fdesc, void *data, size_t n);
void Write( int fdesc, const void *data,size_t n);
void BlockingRead( int fdesc, void *data, size_t n);
pid_t Fork();
void Pipe( int *desc);
void putils_adjust_file_descriptors( pid_t p, int * fds1, int * fds2, int * pfds);
void putils_adjust_file_descriptors_one_way( pid_t p, int * fds, int * pfds, int parent_writes);
size_t tokenize_c_string( char * str, char *delims, char ** tokens, size_t max_ntokens);
char * c_string_to_upper(const char *str);
char * c_string_to_lower( const char *str);
char * join_c_strings(int argc,char **argv);
#ifdef __cplusplus
}
#endif

