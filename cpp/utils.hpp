#ifndef SVM_UTILS_HPP
#define SVM_UTILS_HPP
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cctype>
#include <cstddef>
#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#define TAU 1.e-14

void parse_error(const char *mess);
void quit_error(const char *mess);
FILE *Fopen(const char *name,const char *mode);
size_t Fread(void *ptr,size_t osize,size_t cnt,FILE *fp);
size_t Fwrite(const void *ptr,size_t osize,size_t cnt,FILE *fp);
size_t write_stream(const void *ptr,size_t osize,size_t cnt,std::ostream& os);
size_t read_stream(void *ptr,size_t osize,size_t cnt,std::istream& is);
int explode_string(const std::string& str,const std::string& delims,std::vector<std::string>& tokens);
double dpowi(double x,int i) noexcept;
int parse_bool(const std::string& str) noexcept;
void analyze(int ntp,int ntn,int nfp,int nfn) noexcept;
void write_tdo_file(const char *name,int nvecs,int nfeat,const double *y,const double *vecs);
void write_libsvm_file(const char *name,int nvecs,int nfeat,const double *y,const double *vecs);
void read_tdo_file(const char *name,int *nvecs,int *nfeat,double **y,double **vecs);
void read_libsvm_file(const char *name,int nvecs,int nfeat,double **y,double **vecs);
#endif
