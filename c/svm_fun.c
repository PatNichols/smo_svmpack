#include "svm_fun.h"

inline void analyze(size_t ntp,size_t ntn,size_t nfp,size_t nfn)
{
    const size_t nerrors = nfn + nfp;
    const size_t ncorrect = ntp + ntn;
    const size_t ntotal = ntn + ntp + nfn + nfp;

    const size_t npos = ntp + nfp;
    const size_t nneg = ntn + nfn;
    const size_t nfalse = ntn + nfp;
    const size_t ntrue = ntp + nfn;
    
    double dtp = ntp;
    double dtn = ntn;
    double dfp = nfp;
    double dfn = nfn;
    double dp=npos;
    double dn=nneg;
    double dt=ntrue;
    double df=nfalse;
    double acc = (dtp+dtn)/(dtp+dtn+dfp+dfn);
    double sens = dtp/(dtp+dfn);
    double spec = dtn/(dtn+dfp);
    double prec = dtp/(dtp+dfp);
    double npv = dtn/(dtn+dfn);
    double bal_acc = (sens+spec)*0.5; 
    double f1_score = 2. * sens * prec / ( sens+prec);
    double mcc = (dtp*dtn-dfp*dfn)/sqrt(dt*df*dp*dn);
    fprintf(stderr,"# predictions          = %d\n", ntotal);
    fprintf(stderr,"# true                 = %d\n", ntrue);
    fprintf(stderr,"# false                = %d\n", nfalse);
    fprintf(stderr,"# positive             = %d\n", npos);
    fprintf(stderr,"# negative             = %d\n", nneg);
    fprintf(stderr,"# true-positive        = %d\n", ntp);
    fprintf(stderr,"# true-negative        = %d\n", nfn);
    fprintf(stderr,"# false-positive       = %d\n", nfp);
    fprintf(stderr,"# false-negative       = %d\n", ntn);
    fprintf(stderr,"accuracy               = %le\n", acc);
    fprintf(stderr,"balanced accuracy      = %le\n", bal_acc);
    fprintf(stderr,"# errors               = %d\n", (nfp+nfn));
    fprintf(stderr,"precision              = %le\n", prec);
    fprintf(stderr,"sensitivity(true acc)  = %le\n", sens );
    fprintf(stderr,"specificity(false acc) = %le\n", spec);
    fprintf(stderr,"negative pred. value   = %le\n", npv);        
    fprintf(stderr,"F measure              = %le\n", f1_score );
    fprintf(stderr,"Matthews's Correlation = %le\n", mcc );
}
