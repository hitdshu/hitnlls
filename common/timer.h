#pragma once

#include <chrono>
#include <string>

namespace hitnlls {
namespace common {

class Timer {
public:
    Timer() { is_ticking_ = false; }

    void Tic() { is_ticking_ = true; start_clock_ = std::chrono::high_resolution_clock::now(); }
    double Toc() {
        if (!is_ticking_) {
            return 0;
        }
        is_ticking_ = false;
        end_clock_ = std::chrono::high_resolution_clock::now();
        double t = std::chrono::duration<double, std::milli>(end_clock_ - start_clock_).count();
        return t;
    }

private:
    bool is_ticking_;
    std::chrono::time_point<std::chrono::high_resolution_clock> start_clock_;
    std::chrono::time_point<std::chrono::high_resolution_clock> end_clock_;
};

} // namespace common
} // namespace hitnlls