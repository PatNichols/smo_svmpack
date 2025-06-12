#include "stopwatch.h"

inline stopwatch_t * stopwatch_init() {
    stopwatch_t * s = (stopwatch_t*)Malloc(sizeof(stopwatch_t));
    s->elapsed_time = 0.;
    return s;
}

inline void stopwatch_stop(stopwatch_t *s) {
    struct timespec fn;
    clock_gettime(CLOCK_MONOTONIC,&fn);
    s->elapsed_time += (fn.tv_sec - s->st.tv_sec) + 1.e-9*(fn.tv_nsec-s->st.tv_nsec);
}


