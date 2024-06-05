#include "svmpack_math.hpp"
namespace svmpack {
void analyze(long ntrue_pos, long ntrue_neg, long nfalse_pos, long nfalse_neg)
{
    long n = ntrue_pos + nfalse_pos + ntrue_neg + nfalse_neg;
    long np = ntrue_pos + nfalse_pos;
    long nn = nfalse_neg + ntrue_neg;
    long nt = ntrue_pos + nfalse_neg;
    long nf = ntrue_neg + nfalse_pos;
    double acc = double(ntrue_pos + ntrue_neg)/n;
    double balanced_acc = ( double(ntrue_pos)/np  + double(ntrue_neg)/nn) * 0.5;
    double sensitivity = double(ntrue_pos)/nt;
    double specificity = double(ntrue_neg)/nn;
    double f1 = double(ntrue_pos*2)/double(2*ntrue_pos + nfalse_pos + nfalse_neg);
    double mcc = 0.5 * ( sensitivity + specificity);    
    std::cout << " N = " << n << "\n";
    std::cout << "# errors           = " << (nfalse_pos + nfalse_neg) << "\n";
    std::cout << " accuracy          = " << acc << "\n";
    std::cout << " ntrue_pos         = " << ntrue_pos << "\n";
    std::cout << " ntrue_neg         = " << ntrue_neg << "\n";
    std::cout << " nfalse_pos        = " << nfalse_pos << "\n";
    std::cout << " nfalse_neg        = " << nfalse_neg << "\n";
    std::cout << " sensitivity       = " << sensitivity << "\n";
    std::cout << " specificity       = " << specificity << "\n";
    std::cout << " balanced accuracy = " << balanced_acc << "\n";
    std::cout << " f1 score          = " << f1 << "\n";
    std::cout << " Matthews CC       = " << mcc << "\n";
}
}
