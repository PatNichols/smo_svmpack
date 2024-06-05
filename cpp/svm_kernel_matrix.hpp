#ifndef SVM_KERNEL_MATRIX_HPP
#define SVM_KERNEL_MATRIX_HPP
#include <cmath>
#include <cstdlib>
#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "stopwatch.hpp"
#include "svm_options.hpp"
#include "svm_data.hpp"

namespace svmpack
{
class svm_kernel_matrix
{
private:
    double * cache_rows;
    int64_t * cache_index;
    int64_t cache_size;
    int64_t last;
    double * scale;
    const double * vecs;
    double c1,c2;
    int ktype,kpow;
    int nfeat;
    int kscale;
    int nvecs;

    void precompute_kernel_matrix() noexcept
    {
#pragma omp parallel for
        for (int64_t i=0; i<nvecs; ++i)
        {
            const double * veci = vecs + i * nfeat;
            const double& scali = scale[i];
            double * row = cache_rows + i * nvecs;
            switch (ktype)
            {
            case 0:
                for (int64_t j=0; j<nvecs; ++j)
                {
                    row[j] = scali * scale[j] * (c1*dot(veci,vecs+j*nfeat)+c2);
                }
                break;
            case 1:
                for (int64_t j=0; j<nvecs; ++j)
                {
                    row[j] = scali * scale[j] * std::pow(c1*dot(veci,vecs+j*nfeat)+c2,kpow);
                }
                break;
            case 2:
                for (int64_t j=0; j<nvecs; ++j)
                {
                    row[j] = std::exp(-(c1*dot(veci,vecs+j*nfeat)+c2));
                }
                break;
            case 3:
                for (int64_t j=0; j<nvecs; ++j)
                {
                    row[j] = scali * scale[j] * std::tanh(c1*dot(veci,vecs+j*nfeat)+c2);
                }
                break;
            }
        }
    }

    void init_scale_factors() noexcept
    {
        try
        {
            scale = new double[nvecs];
        }
        catch (...)
        {
            std::cerr << " can not allocate scale factors\n";
            exit(-1);
        }
        if (ktype==2 || !kscale)
        {
#pragma omp parallel for
            for (int i=0; i<nvecs; ++i) scale[i] = 1.0;
        }
        else
        {
            switch (ktype)
            {
            case 0:
#pragma omp parellel for
                for (int i=0; i<nvecs; ++i) scale[i] = c1 * dot(vecs+i*nfeat,vecs+i*nfeat) + c2;
                break;
            case 1:
#pragma omp parellel for
                for (int i=0; i<nvecs; ++i) scale[i] = std::pow(c1 * dot(vecs+i*nfeat,vecs+i*nfeat) + c2,kpow);
                break;
            case 2:
#pragma omp parellel for
                for (int i=0; i<nvecs; ++i) scale[i] = 1.0;
                break;
            case 3:
#pragma omp parellel for
                for (int i=0; i<nvecs; ++i) scale[i] = std::tanh(c1 * dot(vecs+i*nfeat,vecs+i*nfeat) + c2);
                break;
            case 4:
                break;
            }
#pragma omp parellel for
            for (int i=0; i<nvecs; ++i)
            {
                if ( scale[i] > 1.e-15)
                {
                    scale[i] = 1./sqrt(scale[i]);
                }
                else
                {
                    scale[i] = 1.;
                }
            }
        }
    }

    std::size_t find_cache_size(int& precompute_kernel,
                                std::size_t cache_mem,
                                std::size_t cache_size_)
    {
        if (precompute_kernel)
        {
            if ( cache_mem )
            {
                std::size_t nv = (nvecs + 1 ) * sizeof(double);
                std::size_t n_precomp = nvecs * nv;
                if ( n_precomp > cache_mem)
                {
                    precompute_kernel = 0;
                    cache_size = cache_mem/(nv);
                    std::cerr << " cannot precomputed kernel!\n";
                    std::cerr << " cache memory set to " << cache_mem << "\n";
                    std::cerr << " memory needed is    " << n_precomp << "\n";
                    std::cerr << " resetting cache size\n";
                    precompute_kernel = 0;
                }
            }
            else
            {
                cache_size_ = nvecs;
                return cache_size;
            }
        }
        if ( cache_mem )
        {
            std::size_t nv = (nvecs + 1 ) * sizeof(double);
            cache_size_ = cache_mem/nv;
        }
        else
        {
            if ( cache_size == 0)
            {
                cache_size_ = nvecs/2;
            }
        }
        return cache_size_;
    }

    void init_cache(int precompute_kernel,
                    std::size_t cache_mem,
                    std::size_t cache_size_)
    {
        cache_size = find_cache_size(precompute_kernel,cache_mem,cache_size_);
        cache_rows = nullptr;
        cache_index = nullptr;
        if ( precompute_kernel )
        {
            try
            {
                cache_rows = new double[nvecs*nvecs];
                cache_index = nullptr;
                precompute_kernel_matrix();
                return;
            }
            catch (...)
            {
                std::cerr << "cannot precompute kernel matrix : allocation failed\n";
                precompute_kernel = 0;
            }
        }
        if ( cache_mem) {
            std::size_t sz = nvecs*sizeof(double) + sizeof(std::int64_t);
            cache_size = cache_mem/sz;
        }else{
            if ( cache_size_ == 0) {
                cache_size = nvecs/6;
            }else{
                cache_size = cache_size_;
            }
        }
        while ( cache_size >= 2)
        {
            try
            {
                cache_rows = new double[cache_size*nvecs];
                cache_index = new std::int64_t[cache_size];
                break;
            }
            catch (...)
            {
                cache_size >>= 1;
            }
        }
        if ( cache_size < 2)
        {
            // bare mininum for calculation is 2 rows
            try
            {
                cache_rows = new double[2*nvecs];
                cache_index = new std::int64_t[2];
                cache_size = 2;
            }
            catch (...)
            {
                std::cerr << "cannot allocate enough memory for caclulation\n";
                exit(-1);
            }
        }
        std::cerr << "final # of rows to cache = " << cache_size << "\n";
        if ( cache_index )
        {
#pragma omp parallel for
            for (std::int64_t i=0; i<cache_size; ++i) cache_index[i] = -1;
        }
        last=0;
    }

    constexpr double dot(const double *x,const double *y) const noexcept
    {
        double s{0.0};
        for (int i=0; i<nfeat; ++i) s += x[i] * y[i];
        return s;
    }
    constexpr double diff_nrm2(const double *x,const double *y) const noexcept
    {
        double s{0.0};
        for (int i=0; i<nfeat; ++i)
        {
            double t =  x[i] - y[i];
            s += t * t;
        }
        return s;
    }

    const double * calcRow(std::int64_t irow,std::int64_t ivec) noexcept
    {
        const double * veci = vecs + ivec * nfeat;
        const double scali = scale[ivec];
        double * row = cache_rows+irow*nvecs;
        switch ( ktype)
        {
        case 0:
#pragma omp parallel for
            for (int k=0; k<nvecs; ++k)
            {
                row[k] = scale[k] * scali * ( c1 * dot(veci,vecs+k*nfeat) + c2);
            }
            return row;
        case 1:
#pragma omp parallel for
            for (int k=0; k<nvecs; ++k)
            {
                row[k] = scale[k] * scali * std::pow( c1 * dot(veci,vecs+k*nfeat) + c2,kpow);
            }
            return row;
        case 2:
#pragma omp parallel for
            for (int k=0; k<nvecs; ++k)
            {
                row[k] = std::exp( -( c1 * diff_nrm2(veci,vecs+k*nfeat) + c2));
            }
            return row;
        case 3:
#pragma omp parallel for
            for (int k=0; k<nvecs; ++k)
            {
                row[k] = scale[k] * scali * std::tanh( c1 * dot(veci,vecs+k*nfeat) + c2);
            }
            return row;
        }
        return nullptr;
    }
public:

    svm_kernel_matrix(const svm_options& options,
                      const svm_data& data):
        cache_rows(nullptr),
        cache_index(nullptr),
        cache_size(options.cache_size),
        last(0),
        scale(nullptr),
        vecs(data.vectors()),
        c1(options.kc1),
        c2(options.kc2),
        ktype(options.ktype),
        kpow(options.kpow),
        kscale(options.kscale),
        nfeat(data.num_features()),
        nvecs(data.num_vectors())
    {
        if (ktype < 0 || ktype > 3)
        {
            std::cerr << "bad kernel type for svm_kernel_matrix\n";
            std::cerr << "type = " << ktype << "\n";
            exit(-1);
        }
        stopwatch ktimer;
        ktimer.start();
        init_scale_factors();
        init_cache(options.cache_precompute,options.cache_mem,options.cache_size);
        ktimer.stop();
        std::cerr << "kernel matrix initialized\n";
        std::cerr << "time = " << ktimer.time() << "\n";
        std::cerr << "# of cached rows = " << cache_size << "\n";
        std::cerr << "# of features    = " << nfeat << "\n";
        std::cerr << "is scaled " << kscale << "\n";
        std::cerr << "kernel type = " << ktype << "\n";
        std::cerr << " c1 = " << c1 << " c2 = " << c2 << "\n";
        if ( ktype == 1) {
            std::cerr << " kpow = " << kpow << "\n";
        }
    }

    ~svm_kernel_matrix()
    {
        if (cache_index) delete [] cache_index;
        delete [] cache_rows;
        delete [] scale;
    }

    const double * getRow(std::int64_t irow)
    {
        if ( cache_size == nvecs ) return cache_rows + irow * nvecs;
        std::int64_t ifound = -1;
#ifdef _OPENMP
#pragma omp parallel for reduce(max:ifound)
        for (std::int64_t i=0; i<cache_size; ++i)
        {
            if ( cache_index[i] == irow)
            {
                ifound = i;
            }
        }
#else
        for (std::int64_t i=0; i<cache_size; ++i)
        {
            if ( cache_index[i] == irow)
            {
                ifound = i;
                break;
            }
        }
#endif
        if (ifound != -1)
        {
            return cache_rows + ifound * nvecs;
        }
        ifound = last;
        cache_index[last] = irow;
        ++last;
        last %= cache_size;
        return calcRow(ifound,irow);
    }

    constexpr const double * scaleFactors() const noexcept
    {
        return scale;
    }
};
}
#endif