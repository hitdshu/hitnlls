#pragma once

#include "matrix/dense.h"
#include "common/register.h"

namespace hitnlls {

class KernelBase {
public:
    KernelBase() = default;
    virtual ~KernelBase() = default;

    virtual matrix::Vector3d Robustify(double chi) = 0;

    KernelBase(const KernelBase &) = delete;
    KernelBase &operator=(const KernelBase &) = delete;
    KernelBase &operator=(const KernelBase &&) = delete;
};

HITNLLS_REGISTER_REGISTER(KernelBase)
#define HITNLLS_REGISTER_KERNEL(name) \
    HITNLLS_REGISTER_CLASS(KernelBase, name)

class HuberKernel : public KernelBase {
public:
    HuberKernel(double delta = 5) : delta_(delta) {}
    void SetDelta(double delta) { delta_ = delta; }

    virtual matrix::Vector3d Robustify(double chi) override;

protected:
    double delta_;
};

} // namespace hitnlls