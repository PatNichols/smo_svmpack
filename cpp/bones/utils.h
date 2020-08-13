#ifndef SVM_UTILS_H
#define SVM_UTILS_H
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stddef.h>
#include <math.h>
#define MAX_LINE_SIZE 128
#define MAX_TOKENS 128
#define TAU 1.e-14

void parse_error(const char *mess);
void quit_error(const char *mess);
void * Malloc(size_t n);
void * Calloc(size_t n);
void * Realloc(void *ptr,size_t old_size,size_t new_size);
FILE *Fopen(const char *name,const char *mode);
size_t Fread(void *ptr,size_t cnt,size_t osize,FILE *fp);
size_t Fwrite(void *ptr,size_t cnt,size_t osize,FILE *fp);
char ** tokens_init();
void tokens_free(char **tokens);
int64_t find_str(char *str,const char *sub_str);
int64_t find_first_of(char *str,const char *delims,int p);
int64_t find_first_not_of(const char *str,const char *delims,int p);
int explode_string(char *str,const char *delims,char **tokens);
void translate_to_tdo(const char *svm_file,const char *tdo_file);
void translate_to_libsvm(const char *svm_file,const char *tdo_file);
void translate(char *file);
double eval0(const double* v1,const double *v2,int nfeat);
double eval1(const double* v1,const double *v2,int nfeat);
double dpowi(double x,int i);
int parse_bool(char *str);
void analyze(int ntp,int ntn,int nfp,int nfn);

#define MALLOC_PTR(name) ( name *)Malloc(sizeof( name ) )

#endif
