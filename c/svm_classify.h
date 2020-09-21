#ifndef SVM_CLASSIFY_H
#define SVM_CLASSIFY_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#define SVM_USE_OPENMP
#endif
#include "utils.h"
#include "svm_options.h"
#include "svm_data.h"

void svm_classify(svm_options_t *options);

#endif
