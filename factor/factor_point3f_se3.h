#pragma once

#include "factor/base_factor_impl.h"
#include "factor/base_camera.h"

namespace hitnlls {
namespace factor {

class FactorPoint3fSE3 : public BaseFactorImpl<2, 2> {
public:
    virtual ::hitnlls::matrix::Matrix<float, 2, 1> Evaluate() override;
    virtual ::hitnlls::matrix::Matrixxf Jacobian(int nidx) override;

    void SetCamera(BaseCamera *cam) { cam_ = cam; }
protected:
    BaseCamera *cam_;
};

} // namespace factor
} // namespace hitnlls