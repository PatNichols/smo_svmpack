#pragma once
#include <cstdlib>
#include <iostream>
#include <fstream>
#include "svm_data.hpp"
#include "svm_options.hpp"
#include "svm_model.hpp"
namespace svmpack
{
void svm_validate(const svm_options& options);
void svm_classify(const svm_options& options);
}
