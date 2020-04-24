#pragma once

#include "matrix/dense.h"
#include "geometry/se2.h"
#include "type/vertex.h"

namespace hitnlls {

template <typename T> using Matrix33 = matrix::Matrix<T, 3, 3>;

struct VertexSE2 : public Vertex<3, matrix::Matrix33d> {
    virtual void SetZero() override;
    virtual void ApplyPerturbation(const matrix::Vector3d &pert) override;
};

struct VertexSE2Jet : public VertexAudoDiff<3, Matrix33> {
    using BaseType = VertexAudoDiff<3, Matrix33>;
    using PerturbationJetVector = typename BaseType::PerturbationJetVector;

    virtual void SetZero() override {
        estimate_.SetIdentity();
        SetEstimate(estimate_);
    }
    virtual void ApplyPerturbation(const PerturbationJetVector &pert) override {
        geometry::SE2Tangent<geometry::Jetd> left_update(pert[0], pert[1], pert[2]);
        estimate_jet_ = left_update.Exp().ToTransform() * estimate_jet_;
    }
};

} // namespace hitnlls