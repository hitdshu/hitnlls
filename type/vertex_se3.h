#pragma once

#include "matrix/dense.h"
#include "geometry/se3.h"
#include "type/vertex.h"

namespace hitnlls {

template <typename T> using Matrix44 = matrix::Matrix<T, 4, 4>;

struct VertexSE3 : public Vertex<6, matrix::Matrix44d> {
    virtual void SetZero() override;
    virtual void ApplyPerturbation(const matrix::Vector6d &pert) override;
};

struct VertexSE3Jet : public VertexAudoDiff<6, Matrix44> {
    using BaseType = VertexAudoDiff<6, Matrix44>;
    using PerturbationJetVector = typename BaseType::PerturbationJetVector;

    virtual void SetZero() override {
        estimate_.SetIdentity();
        SetEstimate(estimate_);
    }
    virtual void ApplyPerturbation(const PerturbationJetVector &pert) override {
        matrix::Matrix<geometry::Jetd, 3, 1> t = pert.Block(0, 0, 3, 1);
        matrix::Matrix<geometry::Jetd, 3, 1> r = pert.Block(3, 0, 3, 1);
        geometry::SE3Tangent<geometry::Jetd> left_update(t, r);
        estimate_jet_ = left_update.Exp().ToTransform() * estimate_jet_;
    }
};

} // namespace hitnlls