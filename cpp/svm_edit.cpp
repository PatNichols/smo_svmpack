



void svm_edit(
    svm_data& data,
    size_t nfeat_in)
{
    std::size_t nvecs = data_.num_of_vectors();
    std::size_t nfeat_ = data.num_of_features();
    if (nfeat_ == nfeat_in) {
        return;
    }       
//    svmpack::svm_data data_out(nvecs,nfeat_in);
    double * vout = new double[nvecs*nfeat_in];
    const double * vin = data_.vectors();
    if ( nfeat_ < nfeat_in) 
    {
        for (std::size_t k=0;k<nvecs;++k) {
            for (std::size_t j=0;j<nfeat_;++j) 
            {
                vout[k*nfeat_in+j] = vin[k*nfeat_+j];
            }
            for (std::size_t j=nfeat_;j<nfeat_in;++j)
            {
                vout[k*nfeat_in+j] = 0.0;
            }
        }
    }
    else
    {
        for (std::size_t k=0;k<nvecs;++k) {
            for (std::size_t j=0;j<nfeat_in;++j) 
            {
                vout[k*nfeat_in+j] = vin[k*nfeat_+j];
            }
        }
    } 
    const double * tmp = data.vecs;
    data.vecs = vout;
    delete [] tmp;
    data.nfeat = nfeat_in;
}

void svm_edit(
    double * vecs_,
    std::size_t nfeat_old,
    std::size_t nfeat_,
    size_t nvecs)
{

    if (nfeat_ == nfeat_new) {
        return;
    }       
//    svmpack::svm_data data_out(nvecs,nfeat_new);
    double * vout = new double[nvecs*nfeat_new];
    const double * vin = data_.vectors();
    if ( nfeat_ < nfeat_new) 
    {
        for (std::size_t k=0;k<nvecs;++k) {
            for (std::size_t j=0;j<nfeat_;++j) 
            {
                vout[k*nfeat_new+j] = vin[k*nfeat_+j];
            }
            for (std::size_t j=nfeat_;j<nfeat_new;++j)
            {
                vout[k*nfeat_new+j] = 0.0;
            }
        }
    }
    else
    {
        for (std::size_t k=0;k<nvecs;++k) {
            for (std::size_t j=0;j<nfeat_new;++j) 
            {
                vout[k*nfeat_new+j] = vin[k*nfeat_+j];
            }
        }
    } 
    const double * tmp = data.vecs;
    data.vecs = vout;
    delete [] tmp;
}
}
