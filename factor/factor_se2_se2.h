#pragma once

#include "factor/base_factor_impl.h"

namespace hitnlls {
namespace factor {

class FactorSE2SE2 : public BaseFactorImpl<2, 3> {
public:
    virtual ::hitnlls::matrix::Matrix<float, 3, 1> Evaluate() override;
    virtual ::hitnlls::matrix::Matrixxf Jacobian(int nidx) override;
};

} // namespace factor
} // namespace hitnlls