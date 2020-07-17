  
#pragma once

#include <chrono>
#include <string>

namespace nlls {

class Timer {
public:
    explicit Timer();
    void Tic();
    double Toc();
private:
    bool is_ticking_;
    std::chrono::time_point<std::chrono::high_resolution_clock> start_clock_;
    std::chrono::time_point<std::chrono::high_resolution_clock> end_clock_;
};

} // namespace nlls