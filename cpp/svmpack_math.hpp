#pragma once
#include <cstdlib>
#include <iostream>
namespace svmpack
{
constexpr double dot(const double *x, const double *y, int n) noexcept
{
    double sum{0.0};
    for (auto k=0;k<n;++k) sum += x[k] * y[k];
    return sum;
}

constexpr double diff_nrm2(const double *x, const double *y, int n) noexcept
{
    double sum{0.0};
    for (auto k=0;k<n;++k) {
        double t = x[k] - y[k];
        sum += t * t;
    }
    return sum;
}

void analyze(long ntrue_pos, long ntrue_neg, long nfalse_pos, long nfalse_neg);
} // end namespace
