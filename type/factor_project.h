#pragma once

#include "../matrix/dense.h"
#include "../ils/factor.h"
#include "../camera/pinhole.h"

#include "vertex_rn.h"
#include "vertex_se3.h"

namespace nlls {

template <class CameraType>
class FactorProjectOnlyCamera : public Factor<2, CameraType> {
public:
    using ThisType = FactorProjectOnlyCamera;
    using BaseType = Factor<2, CameraType>;
    using JacobianMatrix = typename BaseType::JacobianMatrix;
    explicit FactorProjectOnlyCamera() = default;
    explicit FactorProjectOnlyCamera(const Vector3f &pt3d, const Vector2f &pt2d) : pt3d_(pt3d) { BaseType::SetMeasurement(pt2d); }

    void SetPoint(const Vector3f &pt3d) { pt3d_ = pt3d; }

    virtual void ErrorAndJacobian(bool chi_only = false) override {
        auto &vt = BaseType::GetVertexTuple();
        CameraType *camera = vt.template GetVertex<0>();
        if (chi_only) {
            Vector2f pt_obs = camera->Ray2Pixel(pt3d_);
            BaseType::residual_ = BaseType::sqrt_info_ * (BaseType::measure_ - pt_obs);
            BaseType::chi_ = pow(BaseType::residual_.Norm(), 2);
            return;
        }
        Vector2f pt_obs = camera->Project(pt3d_, nullptr, &BaseType::jacobian_);
        BaseType::residual_ = BaseType::sqrt_info_ * (BaseType::measure_ - pt_obs);
        BaseType::chi_ = pow(BaseType::residual_.Norm(), 2);
        BaseType::jacobians_ = BaseType::sqrt_info_ * BaseType::jacobians_;
    }

private:
    Vector3f pt3d_;
};

class FactorProjectPt3SE3 : public FactorAutoDiff<2, VertexR3Jet, VertexSE3Jet> {
public:
    using ThisType = FactorProjectPt3SE3;
    using BaseType = FactorAutoDiff<2, VertexR3Jet, VertexSE3Jet>;
    using VertexTupleType = typename BaseType::VertexTupleType;
    using EvalVectorTypeJet = typename BaseType::EvalVectorTypeJet;

    void SetCameraParam(const Vector4f &intri);
    virtual EvalVectorTypeJet operator()(VertexTupleType &vars) override final;

private:
    float fx_;
    float fy_;
    float cx_;
    float cy_;
};

} // namespace nlls