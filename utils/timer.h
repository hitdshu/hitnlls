#pragma once

#include <chrono>
#include <string>
#include <iostream>

namespace hitnlls {
namespace utils {

class Timer {
public:
    void Tic() {
        start_ticking_ = true;
        start_ = std::chrono::high_resolution_clock::now();
    }

    double Toc(::std::string desc = "") {
        if(!start_ticking_)
            return 0;
        start_ticking_ = false;
        end_ = ::std::chrono::high_resolution_clock::now();
        double t = ::std::chrono::duration<double, ::std::milli>(end_ - start_).count();
        ::std::cout << "Timer: " << t << " ms in " << desc << ::std::endl;
        return t;
    }

private:
    bool start_ticking_ = false;
    ::std::chrono::time_point<::std::chrono::high_resolution_clock> start_;
    ::std::chrono::time_point<::std::chrono::high_resolution_clock> end_;
};
    
} // namespace utils
} // namespace hitnlls