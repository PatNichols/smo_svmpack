#ifndef SVM_UTILS_H
#define SVM_UTILS_H
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stddef.h>
#include <math.h>
#include <errno.h>

#define MAX_LINE_SIZE 512
#define MAX_TOKENS 128
#define TAU 1.e-14

void * Malloc(size_t n);
void * Calloc(size_t n);
void * Realloc(void *ptr,size_t old_size,size_t new_size);

FILE *Fopen(const char *name,const char *mode);
size_t Fread(void *ptr,size_t size_of_object,size_t cnt,FILE *fp);
size_t Fwrite(void *ptr,size_t size_of_object,size_t cnt,FILE *fp);
ssize_t Getline(char**,size_t*,FILE *fp);

char ** tokens_init();
void tokens_free(char **tokens);

int explode_string(char *str,const char *delims,char **tokens);

double eval0(const double* v1,const double *v2,int nfeat);
double eval1(const double* v1,const double *v2,int nfeat);
double dpowi(double x,int i);

int parse_bool(char *str);

void analyze(int ntp,int ntn,int nfp,int nfn);

void parse_error(const char *mess);
void quit_error(const char *mess);

void write_tdo_file(const char *file_name,int nvecs,int nfeat,double *y,double *vecs);
void write_libsvm_file(const char *file_name,int nvecs,int nfeat,double *y,double *vecs);

void read_tdo_file(const char *file_name,int *nvecs,int *nfeat,double **y,double **vecs);
void read_libsvm_file(const char *file_name,int *nvecs,int *nfeat,double **y,double **vecs);

#define MALLOC_PTR(name) (name*)Malloc( sizeof(name) )

#define MALLOC_ARRAY(name,n) (name*)Malloc( sizeof(name) * n)

#define MALLOC_PTR_PTR(name) (name**)Malloc(sizeof(name*))

#endif
