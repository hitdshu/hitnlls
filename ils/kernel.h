#pragma once

#include "../common/macros.h"
#include "../common/register.h"
#include "../matrix/dense.h"

namespace nlls {

class KernelBase {
public:
    NLLS_NONCOPYABLE(KernelBase)
    KernelBase() = default;
    virtual ~KernelBase() = default;
    virtual Vector2f Robustify(float chi) const = 0;
};

NLLS_REGISTER_REGISTER(KernelBase)
#define NLLS_REGISTER_KERNEL(name) \
    NLLS_REGISTER_REGISTER(KernelBase, name)

class HuberKernel : public KernelBase {
public:
    HuberKernel(float delta = 1) : delta_(delta) {}
    void SetDelta(float delta) { delta_ = delta; }
    virtual Vector2f Robustify(float chi) const override;
private:
    float delta_;
};

class CauchyKernel : public KernelBase {
public:
    CauchyKernel(float c = 1) : c_(c) {}
    void SetC(float c) { c_ = c; }
    virtual Vector2f Robustify(float chi) const override;
private:
    float c_;
};

class GemanMcClureKernel : public KernelBase {
public:
    GemanMcClureKernel(float sigma = 1) : sigma_(sigma) {}
    void SetSigma(float sigma) { sigma_ = sigma; }
    virtual Vector2f Robustify(float chi) const override;
private:
    float sigma_;
};

class ScDcsKernel : public KernelBase {
public:
    ScDcsKernel(float phi = 1) : phi_(phi) {}
    void SetPhi(float phi) { phi_ = phi; }
    virtual Vector2f Robustify(float chi) const override;
private:
    float phi_;
};

} // namespace nlls