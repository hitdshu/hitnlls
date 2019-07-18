#pragma once

#include "factor/base_factor_impl.h"

namespace hitnlls {
namespace factor {

class FactorPoint2fSE2 : public BaseFactorImpl<2, 2> {
public:
    virtual ::hitnlls::matrix::Matrix<float, 2, 1> Evaluate() override;
    virtual ::hitnlls::matrix::Matrixxf Jacobian(int nidx) override;
};

} // namespace factor
} // namespace hitnlls