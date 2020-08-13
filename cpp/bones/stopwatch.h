#ifndef STOPWATCH_H
#define STOPWATCH_H
#include <time.h>
#include <sys/time.h>
#include "utils.h"

typedef struct {
    struct timespec st;
    double elapsed_time;
} stopwatch_t;

stopwatch_t * timer_init();

#define timer_free(s) free((s))
#define timer_start(s) clock_gettime(CLOCK_MONOTONIC,&((s)->st))

void timer_stop(stopwatch_t *s);

#define timer_clear(s) (s)->elapsed_time =0.0

#endif