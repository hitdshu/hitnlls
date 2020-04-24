#pragma once

#include "matrix/dense.h"
#include "type/vertex.h"

namespace hitnlls {

template <typename T> using Vector2 = matrix::Matrix<T, 2, 1>;

struct VertexPt2 : public Vertex<2, Vector2<double>> {
    virtual void SetZero() override;
    virtual void ApplyPerturbation(const matrix::Vector2d &pert) override;
};

struct VertexPt2Jet : public VertexAudoDiff<2, Vector2> {
    using BaseType = VertexAudoDiff<2, Vector2>;
    using PerturbationJetVector = typename BaseType::PerturbationJetVector;

    virtual void SetZero() {
        matrix::Vector2d v;
        v.SetZero();
        SetEstimate(v);
    }
    virtual void ApplyPerturbation(const PerturbationJetVector &pert)  {
        estimate_jet_ += pert;
    }
};

} // namespace hitnlls