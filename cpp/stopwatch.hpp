#ifndef STOPWATCH_HPP
#define STOPWATCH_HPP
#include <ctime>
#include <sys/time.h>

struct stopwatch  {
    struct timespec st;
    double elapsed_time;

    stopwatch():st(),elapsed_time(0.){
    }
    ~stopwatch() {
    }
    void start() noexcept {
        clock_gettime(CLOCK_MONOTONIC,&st);
    }
    void stop() noexcept {
        struct timespec fn;
        clock_gettime(CLOCK_MONOTONIC,&fn);
        elapsed_time += (fn.tv_sec-st.tv_sec) + 1.e-9*(fn.tv_nsec-st.tv_nsec);
    }
    void clear() noexcept { elapsed_time = 0.;}
    double time() const noexcept { return elapsed_time;}
};

#endif
