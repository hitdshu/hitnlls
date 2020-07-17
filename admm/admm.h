#pragma once

#include "question.h"

namespace nlls {
namespace internal {
class Admm {
public:
    explicit Admm(int iter_num = 20) { iter_num_ = iter_num; }

    void SetQuestion(Question *p) { p_ = p; }
    void Solve();

private:
    Question *p_;
    int iter_num_;
};
} // namespace internal
} // namespace nlls