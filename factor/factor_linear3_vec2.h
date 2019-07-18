#pragma once

#include "factor/base_factor_impl.h"
#include "matrix/matrix.h"

namespace hitnlls {
namespace factor {

class FactorLinear3Vec2 : public BaseFactorImpl<1, 3> {
public:
    virtual ::hitnlls::matrix::Matrix<float, 3, 1> Evaluate() override;
    virtual ::hitnlls::matrix::Matrixxf Jacobian(int nidx) override;

    void SetA(const ::hitnlls::matrix::Matrix32f &matA) { matA_ = matA; }
    void Setb(const ::hitnlls::matrix::Vector3f &b) { b_ = b; }

protected:
    ::hitnlls::matrix::Matrix32f matA_;
    ::hitnlls::matrix::Vector3f b_;
};

} // namespace factor
} // namespace hitnlls