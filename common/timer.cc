#include "timer.h"

namespace nlls {

Timer::Timer() { 
    is_ticking_ = false;
}

void Timer::Tic() {
    is_ticking_ = true; 
    start_clock_ = std::chrono::high_resolution_clock::now();
}

double Timer::Toc() {
    if (!is_ticking_) {
        return 0;
    }
    is_ticking_ = false;
    end_clock_ = std::chrono::high_resolution_clock::now();
    double t = std::chrono::duration<double, std::milli>(end_clock_ - start_clock_).count();
    return t;
}

} // namespace nlls