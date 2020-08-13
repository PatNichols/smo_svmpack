#ifdef SVM_USE_OPENMP
#pragma omp parallel shared(max_pair,min_pair,the_max,the_min) private(tid)
    {
        tid = omp_get_thread_num();
        nth = omp_get_num_threads();
        my_max.value= -gmax0;
        my_max.index= -1;
        my_min.value= gmax0;
        my_min.value = -1;
#pragma omp parallel for private(k) shared(nvecs)
        for ( k = 0; k < nvecs; ++k ) {
            ys = smo->y[k] * smo->status[k];
            gx = smo->grad[k];
            if ( ys <= 0. ) {
                if (my_max.value < gx) {
                    my_max.value = gx;
                    my_max.index = k;
                }
            }
            if ( ys >= 0. ) {
                if (my_min.value > gx) {
                    my_min.value = gx;
                    my_min.index = k;
                }
            }
        }
        max_pair[tid].value = my_max.value;
        max_pair[tid].index = my_max.index;
        min_pair[tid].value = my_min.value;
        min_pair[tid].index = my_min.index;
    }
    the_max.value = max_pair[0].value;
    the_max.index = max_pair[0].index;
    for (i=1; i<nth; ++i) {
        if (the_max.value < max_pair[i].value) {
            the_max.value = max_pair[i].value;
            the_max.index = max_pair[i].index;
        }
    }
    the_min.value = min_pair[0].value;
    the_min.index = min_pair[0].index;
    for (i=1; i<nth; ++i) {
        if (the_min.value < min_pair[i].value) {
            the_min.value = min_pair[i].value;
            the_min.index = min_pair[i].index;
        }
    }
#else
