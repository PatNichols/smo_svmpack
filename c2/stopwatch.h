#ifndef STOPWATCH_H
#define STOPWATCH_H
#include <time.h>
#include <sys/time.h>
#include "svm_utils.h"

typedef struct {
    struct timespec st;
    double elapsed_time;
} stopwatch_t;


#define stopwatch_free(s) free((s))
#define stopwatch_start(s) clock_gettime(CLOCK_MONOTONIC,&((s)->st))
#define stopwatch_clear(s) (s)->elapsed_time =0.0
#define stopwatch_time(s) (s)->elapsed_time

stopwatch_t * stopwatch_init();
void stopwatch_stop(stopwatch_t *s);

#endif