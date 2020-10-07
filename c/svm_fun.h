#ifndef SVM_FUN_H
#define SVM_FUN_H
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

inline double svm_dot(const double *v1,const double *v2,int nfeat)
{
    int i;
    double s=0.0;
    for (i=0;i<nfeat;++i) {
        s+=v1[i]*v2[i];
    }
    return s;
}

inline double svm_diff_nrm2(const double *v1,const double *v2,int nfeat)
{
    int i;
    double t;
    double s=0.0;
    for (i=0;i<nfeat;++i) {
        t = v1[i]-v2[i];
        s+=t*t;
    }
    return s;
}

inline double svm_powi(double x,int m)
{
    double w,y,z;
    if (m<0) {
        m = -m;
        x = 1./x;
    }
    switch(m) {
    case 0:
        return 1.0;
    case 1:
        return x;
    case 2:
        return x*x;
    case 3:
        return x*x*x;
    case 4:
        y =x*x;
        return y*y;
    case 5:
        y = x*x;
        return y*y*x;
    case 6:
        y = x*x;
        return y*y*y;
    case 7:
        y = x*x;
        return y*y*y*x;
    case 8:
        y= x*x;
        y= y*y;
        return y*y;
    case 9:
        y = x*x;
        y = y*y;
        return y*y*x;
    case 10:
        y = x*x;
        z = y*y;
        z = y*y;
        return z*y;
    case 11:
        y = x*x;
        z = y*y;
        z = y*y;
        return z*y*x;
    case 12:
        y = x*x;
        z = y*y*y;
        return z*z;
    default:
        y = x*x;
        z = y*y*y;
        z = z*z;
        return z*svm_powi(x,(m-12));
    }
    return 1.e300;
}


void analyze(size_t ntp,size_t ntn,size_t nfp,size_t nfn);

#endif
