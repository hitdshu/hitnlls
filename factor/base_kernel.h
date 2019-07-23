#pragma once

#include "matrix/matrix.h"
#include "common/register.h"

namespace hitnlls {
namespace factor {

class BaseKernel {
public:
    BaseKernel() = default;
    virtual ~BaseKernel() = default;

    virtual ::hitnlls::matrix::Vector3f Robustify(float error) = 0;

    BaseKernel(const BaseKernel &) = delete;
    BaseKernel &operator=(const BaseKernel &) = delete;
};

HITNLLS_REGISTER_REGISTERER(BaseKernel);
#define HITNLLS_REGISTER_KERNEL(name) \
    HITNLLS_REGISTER_CLASS(BaseKernel, name)

} // namespace factor
} // namespace hitnlls