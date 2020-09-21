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
size_t Fread(void *ptr,size_t cnt,size_t osize,FILE *fp);
size_t Fwrite(void *ptr,size_t cnt,size_t osize,FILE *fp);
int explode_string(const std::string& str,const std::string& delims,std::vector<std::string>& tokens);
double dpowi(double x,int i);
int parse_bool(const std::string& str);
void analyze(int ntp,int ntn,int nfp,int nfn);
void write_tdo_file(const char *name,int nvecs,int nfeat,const double *y,const double *vecs);
void write_libsvm_file(const char *name,int nvecs,int nfeat,const double *y,const double *vecs);
#endif
