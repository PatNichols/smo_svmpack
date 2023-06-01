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


template < typename Tp, std::size_t align_val = std::size_t(16) >
struct aligned_allocator {
    typedef Tp value_type;
    typedef Tp* pointer;
    typedef Tp& reference;
    typedef const Tp* const_pointer;
    typedef const Tp& const_reference;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    
    template < class U >
    struct rebind {
        typedef aligned_allocator<U> other;
    };
    
    pointer allocate(size_type n) {
        void * ptr;
#ifdef HAVE_POSIX_MEMALIGN        
        int e = posix_memalign(ptr,align_val,n*sizeof(value_type));
#else
        ptr = malloc(n*sizeof(Tp));
#endif
        if (ptr) return reinterpret_cast< pointer >(ptr);
        std::cerr << "could not allocate " << n << " value_types of size " << sizeof(value_type) << "\n";
        throw std::bad_alloc();
    }
    
    void deallocate(pointer p) {
        free(p);
    }
};


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
