#pragma once

#include "factor/base_kernel.h"

namespace hitnlls {
namespace factor {

class HuberKernel : public BaseKernel {
public:
    HuberKernel(float delta = 5) : delta_(delta) {}

    void SetDelta(float delta) { delta_ = delta; }
    virtual ::hitnlls::matrix::Vector3f Robustify(float error) override;

protected:
    float delta_;
};

} // namespace factor
} // namespace hitnlls